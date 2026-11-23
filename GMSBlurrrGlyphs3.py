# MenuTitle: GMS Blurrr
# -*- coding: utf-8 -*-
__doc__ = """
Blurrr by Gregor Maria Sahl: Gaussian-blurs the outlines of selected glyphs (or only selected paths) and cuts the result at a threshold back into “clean” Bézier paths.
Live preview with sliders for blur, threshold, grunge, grain, smoothness and grid.
Version for Glyphs 3.
"""

import math
import array
import random
import zlib
import traceback

# =====================================================================
#  Geometry core (independent of Glyphs)
# =====================================================================


def marching_squares(rows, w, h, t):
	"""Iso-contours at value t. rows[j][i], j from bottom to top.
	Returns closed polylines in pixel coordinates (sub-pixel accurate)."""
	adj = {}

	def link(a, b):
		adj.setdefault(a, []).append(b)
		adj.setdefault(b, []).append(a)

	for j in range(h - 1):
		r0 = rows[j]
		r1 = rows[j + 1]
		for i in range(w - 1):
			bl = r0[i]
			br = r0[i + 1]
			tr = r1[i + 1]
			tl = r1[i]
			c = (bl >= t) | ((br >= t) << 1) | ((tr >= t) << 2) | ((tl >= t) << 3)
			if c == 0 or c == 15:
				continue
			B = (0, i, j)
			T = (0, i, j + 1)
			L = (1, i, j)
			R = (1, i + 1, j)
			if c == 1 or c == 14:
				link(L, B)
			elif c == 2 or c == 13:
				link(B, R)
			elif c == 3 or c == 12:
				link(L, R)
			elif c == 4 or c == 11:
				link(R, T)
			elif c == 6 or c == 9:
				link(B, T)
			elif c == 7 or c == 8:
				link(L, T)
			else:  # saddle cases 5 and 10: cell average decides
				centerInside = (bl + br + tr + tl) * 0.25 >= t
				if c == 5:
					if centerInside:
						link(L, T); link(B, R)
					else:
						link(L, B); link(R, T)
				else:
					if centerInside:
						link(L, B); link(R, T)
					else:
						link(B, R); link(L, T)

	def edgePoint(e):
		k, i, j = e
		if k == 0:
			a = rows[j][i]; b = rows[j][i + 1]
			return (i + (t - a) / (b - a), float(j))
		a = rows[j][i]; b = rows[j + 1][i]
		return (float(i), j + (t - a) / (b - a))

	loops = []
	visited = set()
	for start in adj:
		if start in visited:
			continue
		visited.add(start)
		loop = [edgePoint(start)]
		prev, cur = None, start
		while True:
			nb = adj[cur]
			if len(nb) != 2:
				break
			nxt = nb[0] if nb[0] != prev else nb[1]
			if nxt == start or nxt in visited:
				break
			visited.add(nxt)
			loop.append(edgePoint(nxt))
			prev, cur = cur, nxt
		if len(loop) >= 4:
			loops.append(loop)
	return loops


# ---------- small vector helpers ----------

def _sub(a, b): return (a[0] - b[0], a[1] - b[1])
def _add(a, b): return (a[0] + b[0], a[1] + b[1])
def _mul(a, s): return (a[0] * s, a[1] * s)
def _dot(a, b): return a[0] * b[0] + a[1] * b[1]
def _cross(a, b): return a[0] * b[1] - a[1] * b[0]
def _len(a): return math.hypot(a[0], a[1])


def _norm(a):
	l = _len(a)
	return (a[0] / l, a[1] / l) if l > 1e-12 else (1.0, 0.0)


def polygon_area(pts):
	s = 0.0
	n = len(pts)
	for k in range(n):
		x0, y0 = pts[k]
		x1, y1 = pts[(k + 1) % n]
		s += x0 * y1 - x1 * y0
	return s * 0.5


def dedupe(pts, eps=1e-4):
	out = []
	for p in pts:
		if not out or _len(_sub(p, out[-1])) > eps:
			out.append(p)
	while len(out) > 1 and _len(_sub(out[0], out[-1])) <= eps:
		out.pop()
	return out


# ---------- curve fitting (Schneider, Graphics Gems) ----------

def _bez(b, t):
	mt = 1.0 - t
	a = mt * mt * mt; c1 = 3 * mt * mt * t; c2 = 3 * mt * t * t; d = t * t * t
	return (a * b[0][0] + c1 * b[1][0] + c2 * b[2][0] + d * b[3][0],
			a * b[0][1] + c1 * b[1][1] + c2 * b[2][1] + d * b[3][1])


def _chord_params(pts):
	u = [0.0]
	for k in range(1, len(pts)):
		u.append(u[-1] + _len(_sub(pts[k], pts[k - 1])))
	tot = u[-1] or 1.0
	return [x / tot for x in u]


def _gen_bezier(pts, u, t1, t2):
	p0, p3 = pts[0], pts[-1]
	c00 = c01 = c11 = x0 = x1 = 0.0
	for p, ui in zip(pts, u):
		mt = 1.0 - ui
		b0 = mt * mt * mt; b1 = 3 * mt * mt * ui; b2 = 3 * mt * ui * ui; b3 = ui * ui * ui
		a1 = _mul(t1, b1); a2 = _mul(t2, b2)
		c00 += _dot(a1, a1); c01 += _dot(a1, a2); c11 += _dot(a2, a2)
		tmp = (p[0] - (p0[0] * (b0 + b1) + p3[0] * (b2 + b3)),
			   p[1] - (p0[1] * (b0 + b1) + p3[1] * (b2 + b3)))
		x0 += _dot(a1, tmp); x1 += _dot(a2, tmp)
	det = c00 * c11 - c01 * c01
	seg = _len(_sub(p3, p0))
	al = ar = 0.0
	if abs(det) > 1e-12:
		al = (x0 * c11 - x1 * c01) / det
		ar = (c00 * x1 - c01 * x0) / det
	eps = 1e-6 * seg
	if al < eps or ar < eps or al > 3 * seg or ar > 3 * seg:
		al = ar = seg / 3.0
	return [p0, _add(p0, _mul(t1, al)), _add(p3, _mul(t2, ar)), p3]


def _newton(b, p, u):
	d = _sub(_bez(b, u), p)
	q1 = [_mul(_sub(b[k + 1], b[k]), 3) for k in range(3)]
	q2 = [_mul(_sub(q1[k + 1], q1[k]), 2) for k in range(2)]
	mt = 1 - u
	d1 = (mt * mt * q1[0][0] + 2 * mt * u * q1[1][0] + u * u * q1[2][0],
		  mt * mt * q1[0][1] + 2 * mt * u * q1[1][1] + u * u * q1[2][1])
	d2 = (mt * q2[0][0] + u * q2[1][0], mt * q2[0][1] + u * q2[1][1])
	num = _dot(d, d1)
	den = _dot(d1, d1) + _dot(d, d2)
	if abs(den) < 1e-12:
		return u
	return min(1.0, max(0.0, u - num / den))


def _max_error(pts, b, u):
	maxd, split = 0.0, len(pts) // 2
	for k in range(1, len(pts) - 1):
		d = _sub(_bez(b, u[k]), pts[k])
		dd = d[0] * d[0] + d[1] * d[1]
		if dd >= maxd:
			maxd, split = dd, k
	return math.sqrt(maxd), split


def fit_cubic(pts, t1, t2, err, maxDepth, depth=0):
	if len(pts) == 2:
		d = _len(_sub(pts[1], pts[0])) / 3.0
		return [[pts[0], _add(pts[0], _mul(t1, d)), _add(pts[1], _mul(t2, d)), pts[1]]]
	u = _chord_params(pts)
	b = _gen_bezier(pts, u, t1, t2)
	e, split = _max_error(pts, b, u)
	if e <= err:
		return [b]
	for _ in range(20 if depth >= maxDepth else 10):
		u = [_newton(b, p, ui) for p, ui in zip(pts, u)]
		b = _gen_bezier(pts, u, t1, t2)
		e, split = _max_error(pts, b, u)
		if e <= err:
			return [b]
	if depth >= maxDepth or len(pts) < 4:
		return [b]
	split = max(1, min(len(pts) - 2, split))
	tc = _norm(_sub(pts[split - 1], pts[split + 1]))
	return (fit_cubic(pts[:split + 1], t1, tc, err, maxDepth, depth + 1) +
			fit_cubic(pts[split:], (-tc[0], -tc[1]), t2, err, maxDepth, depth + 1))


# ---------- extremum detection (with hysteresis against noise) ----------

def _zigzag_extrema(vals, hyst):
	n = len(vals)
	i0 = max(range(n), key=vals.__getitem__)
	res = [i0]
	lookingMin = True
	ext, extI = vals[i0], i0
	for k in range(1, n + 1):
		i = (i0 + k) % n
		v = vals[i]
		if lookingMin:
			if v < ext:
				ext, extI = v, i
			elif v > ext + hyst:
				res.append(extI); lookingMin = False; ext, extI = v, i
		else:
			if v > ext:
				ext, extI = v, i
			elif v < ext - hyst:
				res.append(extI); lookingMin = True; ext, extI = v, i
	return res if len(res) >= 2 else []


def _plateau(vals, e, hf):
	n = len(vals)
	E = vals[e]
	a = b = e
	steps = 0
	while steps < n - 1 and abs(vals[(a - 1) % n] - E) <= hf:
		a = (a - 1) % n; steps += 1
	while steps < n - 1 and abs(vals[(b + 1) % n] - E) <= hf:
		b = (b + 1) % n; steps += 1
	return a, b, steps


def _sign(v):
	return 1.0 if v > 0 else (-1.0 if v < 0 else 0.0)


# ---------- the tracer ----------

def trace_loop(pts, step, tol, maxDepth, cornerDeg, cornerR, lineTol):
	"""Closed polyline -> segments
	('line', p, smooth) / ('curve', c1, c2, p, smooth), starting after the first node."""
	pts = list(pts)
	n = len(pts)
	cum = [0.0]
	for k in range(n):
		cum.append(cum[-1] + _len(_sub(pts[(k + 1) % n], pts[k])))
	total = cum[-1]
	avg = total / n

	def arc(a, b):
		d = cum[b] - cum[a]
		return d if d >= 0 else d + total

	def idist(a, b):
		d = abs(a - b) % n
		return min(d, n - d)

	def between(a, b):  # indices strictly between a and b (forward)
		span = (b - a) % n
		return [(a + j) % n for j in range(1, span)]

	K = max(2, int(round(cornerR / avg)))

	# ---- 1. corners ----
	corners = {}
	Kd = 2 * K  # window for measuring the turning angle (wider than the rounding)
	if cornerDeg < 179.5 and n >= 4 * Kd + 4:
		thr = math.radians(cornerDeg)
		ang = []
		for i in range(n):
			v1 = _norm(_sub(pts[i], pts[(i - Kd) % n]))
			v2 = _norm(_sub(pts[(i + Kd) % n], pts[i]))
			ang.append(math.acos(max(-1.0, min(1.0, _dot(v1, v2)))))
		cand = []
		for i in range(n):
			a = ang[i]
			if a < thr:
				continue
			if all(a > ang[(i - j) % n] for j in range(1, Kd + 1)) and \
					all(a >= ang[(i + j) % n] for j in range(1, Kd + 1)):
				cand.append(i)
		cand.sort()
		m = len(cand)
		for ci, i in enumerate(cand):
			gapP = (i - cand[ci - 1]) % n if m > 1 else n
			gapN = (cand[(ci + 1) % m] - i) % n if m > 1 else n
			mb = max(K + 1, min(3 * K, gapP // 2))
			mf = max(K + 1, min(3 * K, gapN // 2))
			A = pts[(i - K) % n]; A2 = pts[(i - mb) % n]
			B = pts[(i + K) % n]; B2 = pts[(i + mf) % n]
			d1 = _norm(_sub(A, A2))   # direction arriving at the corner
			d2 = _norm(_sub(B2, B))   # direction leaving the corner
			newP = None
			cr = _cross(d1, d2)
			if abs(cr) > 0.08:
				s = _cross(_sub(B, A), d2) / cr
				X = _add(A, _mul(d1, s))
				if _len(_sub(X, pts[i])) < 2.5 * cornerR + step:
					newP = X
			corners[i] = {'tin': d1, 'tout': d2, 'pt': newP}
		for i, c in corners.items():
			if c['pt'] is not None:
				pts[i] = c['pt']

	def nearCorner(i):
		best, bd = None, 2 * K + 1
		for c in corners:
			d = idist(i, c)
			if d < bd:
				best, bd = c, d
		return best

	# ---- 2. extremes and axis-aligned flats ----
	bps = {}
	for i, c in corners.items():
		bps[i] = {'tin': c['tin'], 'tout': c['tout'], 'corner': True, 'prot': True,
				  'lineTo': None, 'skip': K}

	def addSmooth(i, tan, prot=False):
		if i in bps:
			bps[i]['prot'] = bps[i]['prot'] or prot
			return
		bps[i] = {'tin': tan, 'tout': tan, 'corner': False, 'prot': prot,
				  'lineTo': None, 'skip': 0}

	hyst = max(tol * 0.5, step * 0.1)
	hf = max(0.15 * lineTol, 0.1 * tol, 0.02)
	minFlat = max(8.0 * tol, 6.0 * step, 20.0)
	snaps = []  # (indices, axis, value)

	for axis in (0, 1):
		vals = [p[axis] for p in pts]
		for e in _zigzag_extrema(vals, hyst):
			a, b, steps = _plateau(vals, e, hf)
			d = _sub(pts[b], pts[a])
			if steps < n - 2 and arc(a, b) >= minFlat:
				tan = (0.0, _sign(d[1])) if axis == 0 else (_sign(d[0]), 0.0)
				if tan == (0.0, 0.0):
					tan = _norm(d)
				ca, cb = nearCorner(a), nearCorner(b)
				a2 = ca if ca is not None else a
				b2 = cb if cb is not None else b
				if ca is None:
					addSmooth(a2, tan, prot=True)
				if cb is None:
					addSmooth(b2, tan, prot=True)
				bps[a2]['lineTo'] = b2
				E = sum(vals[(a + k) % n] for k in range(steps + 1)) / (steps + 1)
				snaps.append(((a2, b2), axis, E))
			else:
				mid = (a + steps // 2) % n
				if nearCorner(mid) is not None:
					continue
				d = _sub(pts[(mid + 2) % n], pts[(mid - 2) % n])
				tan = (0.0, _sign(d[1])) if axis == 0 else (_sign(d[0]), 0.0)
				if tan == (0.0, 0.0):
					tan = _norm(d)
				addSmooth(mid, tan)

	for (a, b), axis, E in snaps:
		for i in (a, b):
			p = list(pts[i]); p[axis] = E; pts[i] = tuple(p)

	idx = sorted(bps)
	lines = [(i, bps[i]['lineTo']) for i in idx if bps[i]['lineTo'] is not None]

	def inside(i, a, b):
		return 0 < arc(a, i) < arc(a, b)

	idx = [i for i in idx if bps[i]['prot'] or not any(inside(i, a, b) for a, b in lines)]

	# merge nodes that are too close
	mergeDist = max(1.5 * step, tol)
	changed = True
	while changed and len(idx) > 2:
		changed = False
		for k in range(len(idx)):
			a, b = idx[k], idx[(k + 1) % len(idx)]
			if arc(a, b) < mergeDist and bps[a]['lineTo'] != b:
				drop = b if not bps[b]['prot'] else (a if not bps[a]['prot'] else None)
				if drop is not None:
					idx.remove(drop); changed = True
					break

	if len(idx) < 2:
		idx = [0, n // 2]
		bps = {}
		for i in idx:
			t = _norm(_sub(pts[(i + 2) % n], pts[(i - 2) % n]))
			bps[i] = {'tin': t, 'tout': t, 'corner': False, 'prot': False, 'lineTo': None, 'skip': 0}

	m = len(idx)

	def interior(A, B):
		sA = bps[A]['skip']; sB = bps[B]['skip']
		inner = between(A, B)
		if sA + sB >= len(inner):
			return []
		return inner[sA:len(inner) - sB]

	# ---- 3. straight lines ----
	isLine = []
	for k in range(m):
		A, B = idx[k], idx[(k + 1) % m]
		if bps[A]['lineTo'] == B:
			isLine.append(True); continue
		ok = False
		if lineTol > 0 and arc(A, B) >= 3 * step:
			PA, PB = pts[A], pts[B]
			ch = _sub(PB, PA); L = _len(ch)
			if L > 1e-6:
				dev = 0.0
				for j in interior(A, B):
					dev = max(dev, abs(_cross(ch, _sub(pts[j], PA))) / L)
					if dev > lineTol:
						break
				ok = dev <= lineTol
		isLine.append(ok)

	# merge consecutive lines that are really one straight line
	if lineTol > 0:
		changed = True
		while changed and m > 2:
			changed = False
			for k in range(m):
				i = idx[k]
				if bps[i]['corner']:
					continue
				if not (isLine[k - 1] and isLine[k]):
					continue
				A, B = idx[k - 1], idx[(k + 1) % m]
				PA, PB = pts[A], pts[B]
				ch = _sub(PB, PA); L = _len(ch)
				if L < 1e-6:
					continue
				dev = 0.0
				for j in [i] + interior(A, i) + interior(i, B):
					dev = max(dev, abs(_cross(ch, _sub(pts[j], PA))) / L)
				if dev <= lineTol:
					if bps[A]['lineTo'] == i:
						bps[A]['lineTo'] = B
					del idx[k]; del isLine[k]
					m -= 1
					changed = True
					break

	# Two lines meeting at a non-corner node with only a slight bend would give
	# a corner node in the middle of a straight-ish edge. Turn the shorter of
	# the two into a curve instead, so the node becomes a smooth one.
	maxBend = math.radians(30.0)
	for k in range(m):
		i = idx[k]
		if bps[i]['corner'] or not (isLine[k - 1] and isLine[k]):
			continue
		A, B = idx[k - 1], idx[(k + 1) % m]
		d1 = _norm(_sub(pts[i], pts[A]))
		d2 = _norm(_sub(pts[B], pts[i]))
		if math.acos(max(-1.0, min(1.0, _dot(d1, d2)))) > maxBend:
			continue  # a clear kink: keep it as a corner
		if arc(A, i) < arc(i, B):
			isLine[k - 1] = False
		else:
			isLine[k] = False

	# snap almost horizontal / vertical lines
	for k in range(m):
		if not isLine[k]:
			continue
		A, B = idx[k], idx[(k + 1) % m]
		PA, PB = pts[A], pts[B]
		lim = max(lineTol, 0.5)
		if abs(PA[0] - PB[0]) <= lim:
			x = (PA[0] + PB[0]) * 0.5
			pts[A] = (x, PA[1]); pts[B] = (x, PB[1])
		elif abs(PA[1] - PB[1]) <= lim:
			y = (PA[1] + PB[1]) * 0.5
			pts[A] = (PA[0], y); pts[B] = (PB[0], y)

	# smooth nodes next to a line take the line's direction
	smoothFlag = {}
	for k in range(m):
		i = idx[k]
		inLine = isLine[k - 1]
		outLine = isLine[k]
		b = bps[i]
		if b['corner']:
			smoothFlag[i] = False
			continue
		if inLine and outLine:
			smoothFlag[i] = False
			b['tin'] = _norm(_sub(pts[i], pts[idx[k - 1]]))
			b['tout'] = _norm(_sub(pts[idx[(k + 1) % m]], pts[i]))
		elif inLine:
			t = _norm(_sub(pts[i], pts[idx[k - 1]]))
			b['tin'] = b['tout'] = t
			smoothFlag[i] = True
		elif outLine:
			t = _norm(_sub(pts[idx[(k + 1) % m]], pts[i]))
			b['tin'] = b['tout'] = t
			smoothFlag[i] = True
		else:
			smoothFlag[i] = True

	# ---- 4. curves ----
	segs = []
	for k in range(m):
		A, B = idx[k], idx[(k + 1) % m]
		if isLine[k]:
			segs.append(('line', pts[B], smoothFlag[B]))
			continue
		sp = [pts[A]] + [pts[j] for j in interior(A, B)] + [pts[B]]
		t1 = bps[A]['tout']
		t2 = (-bps[B]['tin'][0], -bps[B]['tin'][1])
		pieces = fit_cubic(sp, t1, t2, tol, maxDepth)
		for pi, bz in enumerate(pieces):
			last = pi == len(pieces) - 1
			segs.append(('curve', bz[1], bz[2], bz[3], smoothFlag[B] if last else True))
	return segs


def _normalized(vals):
	n = len(vals)
	if n == 0:
		return vals
	mean = sum(vals) / n
	var = sum((v - mean) * (v - mean) for v in vals) / n
	sd = math.sqrt(var) or 1.0
	return [(v - mean) / sd for v in vals]


def add_grunge(rows, fine, coarse, amount, threshold=0.5):
	"""Adds normalized noise (two octaves) to the field in place.
	fine / coarse: noise fields of the same size; amount: 0..1."""
	if amount <= 0:
		return rows
	w = len(rows[0])
	flat1 = _normalized([v for r in fine for v in r])
	flat2 = _normalized([v for r in coarse for v in r])
	mix = _normalized([a + 0.6 * b for a, b in zip(flat1, flat2)])
	sd = 0.25 * amount
	for j, r in enumerate(rows):
		base = j * w
		for i in range(w):
			f = r[i]
			# grunge only where there is ink nearby: no spatter in empty space
			weight = f / threshold if f < threshold else 1.0
			r[i] = f + weight * sd * mix[base + i]
	return rows


def close_border(rows, value=-1.0):
	"""Forces the outermost ring of the field below any threshold,
	so every contour is closed inside the grid."""
	for i in range(len(rows[0])):
		rows[0][i] = value
		rows[-1][i] = value
	for r in rows:
		r[0] = value
		r[-1] = value
	return rows


def blur_contours(rows, w, h, ox, oy, step, threshold, tol, minSize, cornerDeg, cornerR):
	"""Field -> list of segment lists in font units (specks below minSize removed)."""
	result = []
	minArea = max((2.0 * step) ** 2, minSize * minSize)
	close_border(rows)
	lineTol = max(0.3, 0.5 * tol)
	for loop in marching_squares(rows, w, h, threshold):
		pts = dedupe([(ox + (x + 0.5) * step, oy + (y + 0.5) * step) for x, y in loop])
		if len(pts) < 6 or abs(polygon_area(pts)) < minArea:
			continue
		result.append(trace_loop(pts, step, tol, 14, cornerDeg, cornerR, lineTol))
	return result


# =====================================================================
#  Glyphs part
# =====================================================================

try:
	from GlyphsApp import Glyphs, GSPath, GSNode, LINE, CURVE, OFFCURVE, Message
	IN_GLYPHS = True
except ImportError:
	IN_GLYPHS = False  # running outside Glyphs (tests)

if IN_GLYPHS:
	from AppKit import (NSData, NSApp, NSGraphicsContext, NSColor, NSAffineTransform, NSBezierPath,
						NSNull, NSPoint)
	import Quartz
	try:
		import vanilla
	except ImportError:
		vanilla = None


def _ci(name):
	obj = getattr(Quartz, name, None)
	if obj is None:
		import CoreImage
		obj = getattr(CoreImage, name)
	return obj


_CICONTEXT = [None]


def _ciContext():
	if _CICONTEXT[0] is None:
		opts = {
			_ci("kCIContextWorkingColorSpace"): NSNull.null(),
			_ci("kCIContextWorkingFormat"): _ci("kCIFormatRGBAf"),
		}
		_CICONTEXT[0] = _ci("CIContext").contextWithOptions_(opts)
	return _CICONTEXT[0]


def blur_field(bezierPath, sigma, step):
	"""Rasterizes the path (linear, antialiased), blurs it with Core Image
	and returns the float field: rows (bottom->top), w, h, ox, oy."""
	bounds = bezierPath.bounds()
	pad = 3.0 * sigma + 3.0 * step
	ox = bounds.origin.x - pad
	oy = bounds.origin.y - pad
	w = int(math.ceil((bounds.size.width + 2 * pad) / step)) + 1
	h = int(math.ceil((bounds.size.height + 2 * pad) / step)) + 1

	gray = Quartz.CGColorSpaceCreateWithName(Quartz.kCGColorSpaceLinearGray)
	ctx = Quartz.CGBitmapContextCreate(None, w, h, 8, 0, gray, Quartz.kCGImageAlphaNone)
	Quartz.CGContextSetGrayFillColor(ctx, 0.0, 1.0)
	Quartz.CGContextFillRect(ctx, Quartz.CGRectMake(0, 0, w, h))

	nsctx = NSGraphicsContext.graphicsContextWithCGContext_flipped_(ctx, False)
	NSGraphicsContext.saveGraphicsState()
	NSGraphicsContext.setCurrentContext_(nsctx)
	try:
		tr = NSAffineTransform.transform()
		tr.setTransformStruct_((1.0 / step, 0.0, 0.0, 1.0 / step, -ox / step, -oy / step))
		p = bezierPath.copy()
		p.transformUsingAffineTransform_(tr)
		NSColor.whiteColor().set()
		p.fill()
	finally:
		NSGraphicsContext.restoreGraphicsState()

	image = _ci("CIImage").imageWithCGImage_(Quartz.CGBitmapContextCreateImage(ctx))
	extent = image.extent()
	sigmaPx = sigma / step
	if sigmaPx > 0.05:
		f = _ci("CIFilter").filterWithName_("CIGaussianBlur")
		f.setValue_forKey_(image.imageByClampingToExtent(), "inputImage")
		f.setValue_forKey_(sigmaPx, "inputRadius")
		image = f.valueForKey_("outputImage").imageByCroppingToRect_(extent)

	return _readRows(image, extent, w, h), w, h, ox, oy


def _readRows(image, extent, w, h):
	"""CIImage -> list of rows (bottom->top) with the red channel as floats."""
	linear = Quartz.CGColorSpaceCreateWithName(Quartz.kCGColorSpaceExtendedLinearSRGB)
	cg = _ciContext().createCGImage_fromRect_format_colorSpace_(
		image, extent, _ci("kCIFormatRGBAf"), linear)
	if cg is None:
		raise RuntimeError("Core Image could not render the image.")
	data = bytes(Quartz.CGDataProviderCopyData(Quartz.CGImageGetDataProvider(cg)))
	stride = Quartz.CGImageGetBytesPerRow(cg) // 4
	arr = array.array('f')
	arr.frombytes(data)
	rows = []
	for j in range(h):  # image row 0 = top -> flip
		base = (h - 1 - j) * stride
		rows.append(arr[base:base + 4 * w:4].tolist())  # red channel
	return rows


def _hash_noise(i0, j0, nw, nh, seed):
	"""Deterministic random values 0..255 for lattice cells (i0..i0+nw, j0..j0+nh).
	Each cell's value depends only on its absolute position and the seed."""
	out = bytearray(nw * nh)
	s = (seed * 1442695041) & 0xffffffff
	k = 0
	for jj in range(nh - 1, -1, -1):  # image rows go top -> bottom
		jh = ((j0 + jj) * 668265263 + s) & 0xffffffff
		for ii in range(nw):
			x = ((i0 + ii) * 374761393 + jh) & 0xffffffff
			x = ((x ^ (x >> 13)) * 1274126177) & 0xffffffff
			out[k] = (x ^ (x >> 16)) & 0xff
			k += 1
	return bytes(out)


def noise_field(w, h, ox, oy, step, grain, seed):
	"""Smooth random noise on the same grid as the blur field.
	A lattice of random values with one cell per 'grain' units, anchored to the
	glyph's coordinates, smoothly enlarged to the field's grid. So the pattern
	stays put when Blur or Grid change; grain = size of the bumps."""
	cell = max(grain, step)  # nothing finer than one grid pixel
	i0 = int(math.floor(ox / cell)) - 2
	j0 = int(math.floor(oy / cell)) - 2
	nw = int(math.ceil((w * step) / cell)) + 5
	nh = int(math.ceil((h * step) / cell)) + 5
	data = _hash_noise(i0, j0, nw, nh, seed)
	gray = Quartz.CGColorSpaceCreateDeviceGray()
	provider = Quartz.CGDataProviderCreateWithCFData(NSData.dataWithBytes_length_(data, len(data)))
	cgimg = Quartz.CGImageCreate(nw, nh, 8, 8, nw, gray, Quartz.kCGImageAlphaNone,
								 provider, None, False, Quartz.kCGRenderingIntentDefault)
	img = _ci("CIImage").imageWithCGImage_(cgimg).imageByClampingToExtent()
	# lattice cell (i, j) covers font units [i*cell, (i+1)*cell)
	scale = cell / step
	t = Quartz.CGAffineTransformMake(scale, 0.0, 0.0, scale,
									 (i0 * cell - ox) / step, (j0 * cell - oy) / step)
	img = img.imageByApplyingTransform_(t)
	f = _ci("CIFilter").filterWithName_("CIGaussianBlur")
	f.setValue_forKey_(img, "inputImage")
	f.setValue_forKey_(0.45 * scale, "inputRadius")
	img = f.valueForKey_("outputImage")
	extent = Quartz.CGRectMake(0, 0, w, h)
	return _readRows(img.imageByCroppingToRect_(extent), extent, w, h)


def _makePath(segs):
	path = GSPath()
	smooth = []
	for s in segs:
		if s[0] == 'line':
			path.nodes.append(GSNode(NSPoint(s[1][0], s[1][1]), LINE))
			smooth.append(s[2])
		else:
			path.nodes.append(GSNode(NSPoint(s[1][0], s[1][1]), OFFCURVE))
			path.nodes.append(GSNode(NSPoint(s[2][0], s[2][1]), OFFCURVE))
			path.nodes.append(GSNode(NSPoint(s[3][0], s[3][1]), CURVE))
			smooth += [False, False, s[4]]
	path.closed = True
	for nd, sm in zip(path.nodes, smooth):
		if nd.type != OFFCURVE:
			nd.smooth = bool(sm)
	return path


def _setShapes(layer, shapes):
	for k in range(len(layer.shapes) - 1, -1, -1):
		del layer.shapes[k]
	for s in shapes:
		layer.shapes.append(s)


MAX_PIXELS = 1600


def snapshotSource(layer, selected=None):
	"""Takes the outline to blur BEFORE anything is changed. Components are
	decomposed here, so a composite like Ccedilla is blurred from the original,
	sharp C - even if the C glyph itself gets blurred first.
	Returns (bezierPath, specks)."""
	shapes = []
	if selected:
		shapes = [layer.shapes[k] for k in sorted(selected)]
	else:
		try:
			shapes = list(layer.copyDecomposedLayer().shapes)
		except Exception:
			shapes = []
	bp = NSBezierPath.bezierPath()
	areas = []
	for shape in shapes:
		try:
			bp.appendBezierPath_(shape.bezierPath)
			pts = [(n.position.x, n.position.y) for n in shape.nodes]
		except Exception:
			continue
		if len(pts) >= 3:
			areas.append(abs(polygon_area(pts)))
	if bp.isEmpty() and not selected:
		full = layer.completeBezierPath
		if full is not None:
			bp = full.copy()
	# Specks: everything smaller than MAX_SPECK is removed - but never smaller
	# than half of the smallest original contour, so i-dots and small
	# counters survive.
	specks = MAX_SPECK
	if areas:
		specks = min(MAX_SPECK, math.sqrt(0.5 * min(areas)))
	return bp, specks


def blurLayer(layer, p, seed, selected, cache, source):
	"""Blurs, grunges and re-traces one layer. Returns (paths, nodes).
	selected: indices of the shapes to change (None = whole glyph);
	source: (bezierPath, specks) from snapshotSource()."""
	bp, specks = source
	keep = []
	if selected:
		keep = [shape.copy() for k, shape in enumerate(layer.shapes) if k not in selected]
	if bp is None or bp.isEmpty():
		return 0, 0
	blur, step = p["blur"], p["step"]
	if blur < 0.5 and p["grunge"] <= 0:
		# all effect sliders at zero: leave the outlines exactly as they are
		nodes = sum(len([n for n in layer.shapes[k].nodes if n.type != OFFCURVE])
					for k in range(len(layer.shapes))
					if (not selected or k in selected) and hasattr(layer.shapes[k], "nodes"))
		return len(selected) if selected else len(layer.shapes), nodes
	b = bp.bounds()
	pad = 3.0 * blur + 3.0 * step + (2.0 * p["grain"] if p["grunge"] > 0 else 0.0)
	extent = max(b.size.width, b.size.height) + 2 * pad
	step = max(step, extent / MAX_PIXELS)
	# The rasterized blur and the noise are cached per glyph, so moving
	# Threshold, Grunge or Smoothness doesn't redo the expensive parts.
	if cache is None:
		cache = {}
	fkey = (round(blur, 2), round(step, 3))
	if cache.get("fkey") != fkey:
		cache["field"] = blur_field(bp, blur, step)
		cache["fkey"] = fkey
		cache.pop("nkey", None)
	field, w, h, ox, oy = cache["field"]
	rows = [r[:] for r in field]
	# keep sharp corners as long as the blur is too small to round them visibly
	cornerDeg = 60.0 if blur <= 2.0 * step else 180.0
	cornerR = max(2.0 * step, 1.5 * blur)
	if p["grunge"] > 0:
		nkey = (fkey, round(p["grain"], 1), seed)
		if cache.get("nkey") != nkey:
			cache["noise"] = (noise_field(w, h, ox, oy, step, p["grain"], seed),
							  noise_field(w, h, ox, oy, step, 3.0 * p["grain"], seed + 1))
			cache["nkey"] = nkey
		fine, coarse = cache["noise"]
		add_grunge(rows, fine, coarse, p["grunge"], p["threshold"])
	paths = [_makePath(segs) for segs in
			 blur_contours(rows, w, h, ox, oy, step, p["threshold"], p["tol"], specks,
						   cornerDeg, cornerR)]
	_setShapes(layer, keep + paths)
	layer.correctPathDirection()
	layer.roundCoordinates()
	nodes = sum(len([n for n in pa.nodes if n.type != OFFCURVE]) for pa in paths)
	return len(paths), nodes


MAX_SPECK = 100.0  # islands / holes smaller than this (units, square side) are removed

PREF = "com.claude.blurrr."
#        key          title              lo     hi     default  left hint             right hint
SLIDERS = (
	("blur",      "Blur (units)",     0.0, 150.0, 20.0, "sharp",               "soft"),
	("threshold", "Threshold",        0.05,  0.95, 0.5, "bold",                "light"),
	("grunge",    "Strength",         0.01,  1.0,  0.01, "clean",              "grungy"),
	("grain",     "Grain (units)",    1.0,  80.0,  6.0, "fine",                "rough"),
	("tol",       "Smoothness",       0.2,   8.0,  1.0, "exact",               "smooth"),
	("step",      "Grid (units/px)",  1.0,  10.0,  3.0, "hi-res",              "lo-res"),
)
SECTIONS = (("Shape", ("blur", "threshold")),
			("Grunge", ("grunge", "grain")),
			("Output", ("tol", "step")))


class BlurDialog(object):

	def __init__(self):
		self.font = Glyphs.font
		if self.font is None:
			Message("No font open.", "GMS Blurrr")
			return
		self.layers = [l for l in self.font.selectedLayers if l.parent is not None]
		if not self.layers:
			Message("Please select some glyphs.", "GMS Blurrr")
			return
		self.originals = [[s.copy() for s in l.shapes] for l in self.layers]
		# selected paths per layer (None = no selection -> whole glyph)
		self.selection = []
		for l in self.layers:
			sel = set()
			for k, shape in enumerate(l.shapes):
				try:
					if any(n.selected for n in shape.nodes):
						sel.add(k)
				except Exception:
					pass  # components etc.
			self.selection.append(sel or None)
		# outlines to blur, taken now - before any glyph is changed
		self.sources = [snapshotSource(l, sel) for l, sel in zip(self.layers, self.selection)]
		self.applied = False
		self.caches = {}

		for sdef in SLIDERS:
			Glyphs.registerDefault(PREF + sdef[0], sdef[4])
		for k, v in (("seed", 0),):
			Glyphs.registerDefault(PREF + k, v)
		self.seed = int(Glyphs.defaults[PREF + "seed"])

		W = 380
		H = 12 + len(SECTIONS) * 20 + (len(SECTIONS) - 1) * 10 + len(SLIDERS) * 40 + 26 + 18 + 34
		self.w = vanilla.FloatingWindow((W, H), "GMS Blurrr")
		self.sliders = {}
		self.hints = {}
		defs = {sdef[0]: sdef for sdef in SLIDERS}
		y = 12
		for si, (section, keys) in enumerate(SECTIONS):
			if si:
				setattr(self.w, "line%d" % si, vanilla.HorizontalLine((15, y, -15, 1)))
				y += 10
			setattr(self.w, "head%d" % si, vanilla.TextBox((15, y - 1, 60, 17), section, sizeStyle="small"))
			if section == "Grunge":
				# off every time the window opens
				self.w.useGrunge = vanilla.CheckBox((66, y - 3, 22, 20), "", value=False,
													callback=self.grungeToggled, sizeStyle="small")
			y += 20
			for key in keys:
				_, title, lo, hi, _, hintLo, hintHi = defs[key]
				setattr(self.w, key + "Title", vanilla.TextBox((15, y + 1, 115, 17), title, sizeStyle="small"))
				sl = vanilla.Slider((130, y, -60, 19), minValue=lo, maxValue=hi,
									value=self._prefValue(key), callback=self.changed, sizeStyle="small")
				setattr(self.w, key + "Slider", sl)
				setattr(self.w, key + "Value", vanilla.TextBox((-55, y + 1, -15, 17), "", sizeStyle="small", alignment="right"))
				lo_ = vanilla.TextBox((132, y + 19, 110, 14), hintLo, sizeStyle="mini")
				hi_ = vanilla.TextBox((-180, y + 19, -60, 14), hintHi, sizeStyle="mini", alignment="right")
				setattr(self.w, key + "HintLo", lo_)
				setattr(self.w, key + "HintHi", hi_)
				self.sliders[key] = sl
				self.hints[key] = (lo_, hi_)
				y += 40
			if section == "Grunge":
				self.w.shuffle = vanilla.Button((-115, y - 2, -15, 20), "Shuffle",
												callback=self.shuffleCallback, sizeStyle="small")
				y += 26
		y += 4
		self.w.lineEnd = vanilla.HorizontalLine((15, y, -15, 1))
		y += 14
		self.w.cancel = vanilla.Button((-225, y, -120, 20), "Cancel", callback=self.cancelCallback, sizeStyle="small")
		self.w.apply = vanilla.Button((-115, y, -15, 20), "Apply", callback=self.applyCallback, sizeStyle="small")
		self.w.setDefaultButton(self.w.apply)
		self.w.status = vanilla.TextBox((15, y + 3, -250, 14), "", sizeStyle="mini")
		self.w.bind("close", self.windowClosed)
		self.updateLabels()

		# The panel must not take keyboard focus, so the space bar keeps
		# showing the black preview fill in the edit view.
		for ctrl in list(self.sliders.values()) + [self.w.useGrunge, self.w.shuffle, self.w.cancel, self.w.apply]:
			try:
				ctrl.getNSView().setRefusesFirstResponder_(True)
			except Exception:
				pass
		previousKey = NSApp().keyWindow()
		self.w.open()
		try:
			self.w.getNSWindow().setBecomesKeyOnlyIfNeeded_(True)
		except Exception:
			pass
		if previousKey is not None:
			previousKey.makeKeyWindow()
		self.preview()

	# ---------- parameters ----------

	def _prefValue(self, key):
		v = float(Glyphs.defaults[PREF + key])
		if key == "threshold" and v > 1.0:  # a test version stored 0 - 100
			v = 0.05 + 0.9 * v / 100.0
		return v

	def params(self):
		p = {k: sl.get() for k, sl in self.sliders.items()}
		if not self.w.useGrunge.get():
			p["grunge"] = 0.0
		return p

	def updateLabels(self):
		p = self.params()
		self.w.blurValue.set("%.0f" % p["blur"])
		self.w.thresholdValue.set("%.2f" % p["threshold"])
		on = bool(self.w.useGrunge.get())
		self.w.grungeValue.set("%d%%" % max(1, round(self.sliders["grunge"].get() * 100)))
		# Strength, Grain and Shuffle only work with "Use grunge" ticked
		self.sliders["grunge"].enable(on)
		self.sliders["grain"].enable(on)
		self.w.shuffle.enable(on)
		for name in ("grungeTitle", "grungeValue", "grungeHintLo", "grungeHintHi",
					 "grainTitle", "grainValue", "grainHintLo", "grainHintHi"):
			try:
				getattr(self.w, name).getNSTextField().setTextColor_(
					NSColor.labelColor() if on else NSColor.disabledControlTextColor())
			except Exception:
				pass
		self.w.grainValue.set("%.0f" % p["grain"])
		self.w.tolValue.set("%.1f" % p["tol"])
		self.w.stepValue.set("%.1f" % p["step"])

	def savePrefs(self):
		for sdef in SLIDERS:
			Glyphs.defaults[PREF + sdef[0]] = float(self.sliders[sdef[0]].get())
		Glyphs.defaults[PREF + "seed"] = self.seed

	def glyphSeed(self, layer):
		# every glyph gets its own pattern; Shuffle changes all of them
		return (zlib.crc32(layer.parent.name.encode("utf-8")) + 7919 * self.seed) & 0x7fffffff

	# ---------- workflow ----------

	def restore(self):
		try:
			self.w.status.set("")
		except Exception:
			pass
		for layer, orig in zip(self.layers, self.originals):
			um = self._undoManager(layer)
			if um: um.disableUndoRegistration()
			try:
				_setShapes(layer, [s.copy() for s in orig])
			finally:
				if um: um.enableUndoRegistration()
		self._redraw()

	def _undoManager(self, layer):
		try:
			return layer.parent.undoManager()
		except Exception:
			return None

	def run(self, withUndo):
		p = self.params()
		if p["grunge"] <= 0.005:
			p["grunge"] = 0.0
		nPaths = nNodes = 0
		errors = []
		self.font.disableUpdateInterface()
		try:
			for layer, orig, sel, src in zip(self.layers, self.originals, self.selection, self.sources):
				um = self._undoManager(layer)
				if um: um.disableUndoRegistration()
				try:
					_setShapes(layer, [s.copy() for s in orig])
				finally:
					if um: um.enableUndoRegistration()
				if withUndo:
					layer.parent.beginUndo()
				elif um:
					um.disableUndoRegistration()
				try:
					a, b = blurLayer(layer, p, self.glyphSeed(layer), sel,
									 self.caches.setdefault(id(layer), {}), src)
					nPaths += a; nNodes += b
				except Exception:
					errors.append(layer.parent.name)
					print("Blurrr - error in /%s:" % layer.parent.name)
					print(traceback.format_exc())
				finally:
					if withUndo:
						layer.parent.endUndo()
					elif um:
						um.enableUndoRegistration()
		finally:
			self.font.enableUpdateInterface()
		self._redraw()
		if errors:
			self.w.status.set("Error in /%s - details in the Macro window" % ", /".join(errors))
			try:
				Glyphs.showMacroWindow()
			except Exception:
				pass
		else:
			note = " (selected paths only)" if any(self.selection) else ""
			self.w.status.set("%d paths, %d nodes%s" % (nPaths, nNodes, note))

	def _redraw(self):
		tab = self.font.currentTab
		if tab is None:
			return
		for name in ("redraw", "forceRedraw"):  # Glyphs 4 / Glyphs 3
			fn = getattr(tab, name, None)
			if fn is not None:
				try:
					fn()
				except Exception:
					pass
				return

	def preview(self):
		self.run(withUndo=False)

	def changed(self, sender):
		self.updateLabels()
		# while dragging only update the labels; compute on mouse-up
		ev = NSApp().currentEvent()
		if ev is not None and ev.type() == 6:  # NSEventTypeLeftMouseDragged
			return
		self.preview()

	def grungeToggled(self, sender):
		# Strength keeps its last used value (0 if never used)
		self.updateLabels()
		self.preview()

	def shuffleCallback(self, sender):
		self.seed += 1
		self.preview()

	def cancelCallback(self, sender):
		# same as the close button: restores the original outlines
		self.w.close()

	def applyCallback(self, sender):
		self.savePrefs()
		self.run(withUndo=True)
		self.applied = True
		self.w.close()

	def windowClosed(self, sender):
		if not self.applied:
			self.savePrefs()
			self.restore()


if IN_GLYPHS:
	Glyphs.clearLog()
	if vanilla is None:
		Message("Blurrr needs the Vanilla module.\n\n"
				"Install it in Window > Plugin Manager > Modules (also Python, if it isn’t installed yet) "
				"and restart Glyphs.", "GMS Blurrr")
	else:
		try:
			BlurDialog()
		except Exception:
			print("GMS Blurrr could not open:")
			print(traceback.format_exc())
			try:
				Glyphs.showMacroWindow()
			except Exception:
				pass
			Message("Blurrr could not open. Details are in the Macro window.", "GMS Blurrr")
