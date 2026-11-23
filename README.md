# Blurrr

A Python script for the [Glyphs](https://glyphsapp.com) font editor that blurs glyph outlines.
Blur the shapes, set a threshold, add some grunge – and get soft, melted or rough versions of your letters right inside Glyphs.

By Gregor Maria Sahl

---

## Download

**[⬇ Download Blurrr (ZIP)](https://github.com/gregormaria/Blurrr/archive/refs/heads/main.zip)**

The ZIP contains two versions:

| File | For |
|---|---|
| `GMSBlurrrGlyphs3.py` | Glyphs 3 |
| `GMSBlurrrGlyphs4.py` | Glyphs 4 |

Only install the file that matches your Glyphs version.

## Installation

1. Unzip the download.
2. In Glyphs, choose **Script → Open Scripts Folder**.
3. Copy the matching `.py` file into that folder.
4. Hold **Option (⌥)** and choose **Script → Reload Scripts** (or restart Glyphs).
5. Blurrr now appears in the **Script** menu.

## How to use

1. Select one or more glyphs.
2. Run **Script → GMS Blurrr**. The glyphs update live as you move the sliders.
3. **Shape**
   - **Blur** – how strongly the outlines are softened (sharp ↔ soft)
   - **Threshold** – where the blurred shape is cut back into an outline (bold ↔ light)
4. **Grunge** – tick the checkbox to add rough, irregular edges
   - **Strength** – how much texture is added (clean ↔ grungy)
   - **Grain** – size of the texture (fine ↔ rough)
   - **Shuffle** – try a new random variation
5. **Output** – controls how the result is turned back into clean Bézier paths
   - **Smoothness** – exact (follows every detail, more nodes) ↔ smooth (fewer nodes, softer curves)
   - **Grid** – resolution of the calculation: hi-res (more detail, slower) ↔ lo-res (faster, coarser)
   
   The status line at the bottom shows how many paths and nodes the result has.
6. Click **Apply** to keep the result, or **Cancel** to restore the original outlines.

Tip: Work on a copy of your glyphs or a duplicate layer, so you can always go back to the original outlines.

## License

Released under the MIT License – free to use, modify and share. See [LICENSE](LICENSE).

## Feedback

Found a bug or have an idea? Open an [issue](https://github.com/gregormaria/Blurrr/issues) here on GitHub.
