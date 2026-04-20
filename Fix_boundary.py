"""
TIF Column/Row Crop-and-Shift Tool
===================================
Overview
--------
This script removes a contiguous gap of columns and a contiguous gap of rows
from a single-band GeoTIFF, then shifts the right/bottom portion of the raster
inward so that the two halves are seamlessly joined.

Concretely, the two edits are:

  Column edit
    Source columns  0 .. 43454          -> kept as-is  (left block)
    Source columns  43455 .. 48045      -> DISCARDED    (4,591 columns removed)
    Source columns  48046 .. max_col    -> shifted left by 4,591 pixels

  Row edit
    Source rows     0 .. 39478          -> kept as-is  (top block)
    Source rows     39479 .. 42163      -> DISCARDED    (2,685 rows removed)
    Source rows     42164 .. max_row    -> shifted up   by 2,685 pixels

Output dimensions
    new_width  = src_width  - 4,591
    new_height = src_height - 2,685

Inherited attributes (all preserved in the output file)
    * Single-band pixel data
    * ColorMap / palette          (via rasterio  write_colormap)
    * Raster Attribute Table      (via GDAL       SetDefaultRAT)
    * CRS, GeoTransform, nodata, dtype, compression, and all other
      profile fields copied verbatim from the source
    * Dataset-level and band-level metadata tags

Memory strategy
    Pixel data is processed in horizontal strips of `chunk_rows` rows so that
    arbitrarily large files (tens of GB) can be handled without loading the
    entire raster into RAM.  The default strip height is 512 rows; reduce it
    with --chunk-rows if memory is tight.

Usage
-----
    python tif_shift.py  input.tif  output.tif  [--chunk-rows 512]

Dependencies
------------
    pip install rasterio numpy
    # For RAT copying (optional):
    conda install -c conda-forge gdal
"""

import argparse
import math

import numpy as np
import rasterio
from rasterio.transform import Affine

# GDAL is only needed for writing the Raster Attribute Table (RAT).
# rasterio does not expose a RAT write API, so we fall back to osgeo.gdal.
# If gdal is not installed the script still runs; RAT copying is simply skipped.
try:
    from osgeo import gdal
    gdal.UseExceptions()
    HAS_GDAL = True
except ImportError:
    HAS_GDAL = False


# ---------------------------------------------------------------------------
# Edit parameters
# All pixel indices are 0-based and follow rasterio convention:
#   column 0 = leftmost pixel,  row 0 = topmost pixel.
# ---------------------------------------------------------------------------

COL_SRC_START = 48046   # First source column of the right block to be shifted
COL_DST_START = 43455   # Destination column where the right block will land
COL_GAP       = COL_SRC_START - COL_DST_START   # Columns removed = 4,591

ROW_SRC_START = 42164   # First source row of the bottom block to be shifted
ROW_DST_START = 39479   # Destination row where the bottom block will land
ROW_GAP       = ROW_SRC_START - ROW_DST_START   # Rows removed = 2,685


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def build_output_transform(src_transform: Affine, src_width: int, src_height: int):
    """
    Compute the output raster dimensions and geo-transform.

    Because the left block (columns 0 .. COL_DST_START-1) is unchanged,
    the top-left corner of the output raster coincides with the top-left
    corner of the source raster.  The geo-transform can therefore be
    reused verbatim; only the pixel dimensions change.

    Parameters
    ----------
    src_transform : Affine
        Affine geo-transform of the source raster.
    src_width, src_height : int
        Pixel dimensions of the source raster.

    Returns
    -------
    out_transform : Affine   (identical to src_transform)
    out_width     : int
    out_height    : int
    """
    out_width  = src_width  - COL_GAP
    out_height = src_height - ROW_GAP
    return src_transform, out_width, out_height


def src_cols_for_out_block(out_col_start: int, out_col_count: int, src_width: int):
    """
    Map a range of output columns to the corresponding source column segments.

    The output column space is divided into two regions:

        Region A  out_col < COL_DST_START
                  Maps 1-to-1 to the same source column (left block, unchanged).

        Region B  out_col >= COL_DST_START
                  Maps to source column = out_col + COL_GAP (right block, shifted).

    Parameters
    ----------
    out_col_start : int
        First output column of the strip being processed.
    out_col_count : int
        Number of output columns in the strip.
    src_width : int
        Total width of the source raster (used as a boundary guard).

    Returns
    -------
    list of (src_col_start, count, out_buf_offset)
        Each tuple describes one contiguous read from the source:
          src_col_start  -- first source column to read
          count          -- number of columns
          out_buf_offset -- offset into the output buffer where this segment lands
    """
    segments = []
    out_end = out_col_start + out_col_count

    # --- Region A: left block (identity mapping) ---
    a_out_start = out_col_start
    a_out_end   = min(out_end, COL_DST_START)
    if a_out_start < a_out_end:
        segments.append((
            a_out_start,                    # src_col_start (== out_col)
            a_out_end - a_out_start,        # count
            a_out_start - out_col_start,    # offset in output buffer
        ))

    # --- Region B: right block (shifted left by COL_GAP) ---
    b_out_start = max(out_col_start, COL_DST_START)
    b_out_end   = out_end
    if b_out_start < b_out_end:
        src_b_start = b_out_start + COL_GAP
        src_b_end   = b_out_end   + COL_GAP
        if src_b_start < src_width:         # guard against out-of-bounds
            src_b_end = min(src_b_end, src_width)
            segments.append((
                src_b_start,
                src_b_end - src_b_start,
                b_out_start - out_col_start,
            ))

    return segments


def out_row_segments(src_row_start: int, src_row_count: int, src_height: int):
    """
    For a source row strip, determine which rows are kept and where they land
    in the output raster.

    Source rows fall into three zones:

        Zone 1  src_row < ROW_DST_START
                Kept; output row = source row  (top block, unchanged).

        Zone 2  ROW_DST_START <= src_row < ROW_SRC_START
                Discarded -- this is the gap that is removed.

        Zone 3  src_row >= ROW_SRC_START
                Kept; output row = source row - ROW_GAP  (bottom block, shifted up).

    Parameters
    ----------
    src_row_start : int
        First source row of the current processing chunk.
    src_row_count : int
        Number of rows in the current chunk.
    src_height : int
        Total height of the source raster (kept for API symmetry).

    Returns
    -------
    list of (src_buf_offset, out_row, n_rows)
        src_buf_offset -- row offset within the chunk's source data array
        out_row        -- first destination row in the output file
        n_rows         -- number of rows in this segment
    """
    segments = []
    src_end = src_row_start + src_row_count

    # --- Zone 1: top block (identity mapping) ---
    z1_start = src_row_start
    z1_end   = min(src_end, ROW_DST_START)
    if z1_start < z1_end:
        segments.append((
            z1_start - src_row_start,   # offset into chunk array
            z1_start,                   # output row index
            z1_end - z1_start,          # count
        ))

    # --- Zone 2: discarded gap -- no entry added ---

    # --- Zone 3: bottom block (shifted up by ROW_GAP) ---
    z3_start = max(src_row_start, ROW_SRC_START)
    z3_end   = src_end
    if z3_start < z3_end:
        segments.append((
            z3_start - src_row_start,
            z3_start - ROW_GAP,         # output row = src_row - ROW_GAP
            z3_end - z3_start,
        ))

    return segments


# ---------------------------------------------------------------------------
# Attribute helpers
# ---------------------------------------------------------------------------

def read_colormap_safe(src: rasterio.DatasetReader, band: int = 1):
    """
    Read the ColorMap from the given band.

    Returns a dict {pixel_value: (R, G, B, A)} if a colormap exists,
    or None if the band carries no colormap.  rasterio raises an exception
    when no colormap is present, which is silently caught here.
    """
    try:
        return src.colormap(band)
    except Exception:
        return None


def copy_rat(src_path: str, dst_path: str, band: int = 1) -> bool:
    """
    Copy the Raster Attribute Table (RAT) from the source file to the
    destination file using GDAL.

    rasterio does not provide a write API for RATs, so this function opens
    the already-written output file in update mode via osgeo.gdal and calls
    SetDefaultRAT().  The destination file must already be closed by rasterio
    before this function is invoked.

    Parameters
    ----------
    src_path : str   Path to the source GeoTIFF.
    dst_path : str   Path to the destination GeoTIFF (must be closed).
    band     : int   1-based band index (default 1).

    Returns
    -------
    True  if the RAT was successfully copied.
    False if the source has no RAT, GDAL is unavailable, or an error occurred.
    """
    if not HAS_GDAL:
        print("  [WARNING] osgeo.gdal not found -- RAT copy skipped.\n"
              "            To enable: conda install -c conda-forge gdal")
        return False

    src_ds = gdal.Open(src_path, gdal.GA_ReadOnly)
    dst_ds = gdal.Open(dst_path, gdal.GA_Update)
    try:
        rat = src_ds.GetRasterBand(band).GetDefaultRAT()
        if rat is None or rat.GetRowCount() == 0:
            print("  [INFO] Source file has no RAT -- skipping.")
            return False
        dst_ds.GetRasterBand(band).SetDefaultRAT(rat)
        dst_ds.FlushCache()
        print(f"  RAT copied successfully "
              f"({rat.GetRowCount()} rows x {rat.GetColumnCount()} columns)")
        return True
    finally:
        # Setting to None triggers GDAL's reference-count destructor,
        # which flushes pending writes and closes the file handles.
        src_ds = None
        dst_ds = None


# ---------------------------------------------------------------------------
# Main processing routine
# ---------------------------------------------------------------------------

def process(input_path: str, output_path: str, chunk_rows: int = 512) -> None:
    """
    Read the source GeoTIFF, apply the column and row crop-and-shift edits,
    and write the result to output_path while preserving all metadata.

    Processing pipeline
    -------------------
    1. Open source; read dimensions, dtype, profile, colormap, and tags.
    2. Create the output file with an updated profile (new width / height).
    3. Iterate over the source in horizontal strips of `chunk_rows` rows.
       For each strip:
         a. Read the full-width strip from the source.
         b. Determine which rows are kept and their output positions
            (via out_row_segments).
         c. For each kept row sub-strip, build an output-width zero buffer
            and fill it from the two column segments (via src_cols_for_out_block).
         d. Write the buffer to the correct window in the output file.
    4. Write the ColorMap before closing the output file (rasterio requirement:
       write_colormap must be called while the dataset is still open).
    5. After the output file is closed, copy the RAT via GDAL in update mode.

    Parameters
    ----------
    input_path  : str  Path to the source single-band GeoTIFF.
    output_path : str  Path where the processed GeoTIFF will be written.
    chunk_rows  : int  Number of source rows processed per iteration.
                       Larger values are faster but consume more RAM.
                       Default: 512 rows ~ 512 x src_width x bytes_per_pixel.
    """
    print(f"Source file : {input_path}")
    with rasterio.open(input_path) as src:
        src_width  = src.width
        src_height = src.height
        band_count = src.count
        dtype      = src.dtypes[0]
        transform  = src.transform

        print(f"  Dimensions : {src_width} cols x {src_height} rows, "
              f"{band_count} band(s), dtype={dtype}")
        assert band_count == 1, (
            f"This script targets single-band TIFFs; source has {band_count} bands."
        )

        # Validate that the specified pixel indices lie within the source raster.
        assert COL_DST_START < COL_SRC_START <= src_width, (
            f"Column index out of range: COL_SRC_START={COL_SRC_START}, "
            f"src_width={src_width}"
        )
        assert ROW_DST_START < ROW_SRC_START <= src_height, (
            f"Row index out of range: ROW_SRC_START={ROW_SRC_START}, "
            f"src_height={src_height}"
        )

        out_transform, out_width, out_height = build_output_transform(
            transform, src_width, src_height
        )

        print(f"  Output     : {out_width} cols x {out_height} rows")
        print(f"  Col gap    : [{COL_DST_START}, {COL_SRC_START})  ->  {COL_GAP} columns removed")
        print(f"  Row gap    : [{ROW_DST_START}, {ROW_SRC_START})  ->  {ROW_GAP} rows removed")

        # ------------------------------------------------------------------
        # Collect attributes to be inherited by the output file
        # ------------------------------------------------------------------

        # ColorMap maps integer pixel values to (R, G, B, A) tuples.
        # Typical in classified rasters (land cover, geology, soil type, etc.).
        colormap = read_colormap_safe(src, band=1)
        print(f"  ColorMap   : "
              f"{'found, ' + str(len(colormap)) + ' entries' if colormap else 'none'}")

        # Dataset-level tags include GDAL standard keys (e.g. AREA_OR_POINT)
        # as well as any custom key-value metadata from the producer.
        dataset_tags = src.tags()

        # Band-level tags may carry per-band descriptions, units, scale/offset, etc.
        band_tags = src.tags(1)

        # ------------------------------------------------------------------
        # Build the output profile
        # ------------------------------------------------------------------
        # src.profile already contains driver, CRS, nodata, dtype, and the
        # original compression settings.  We override only the fields that
        # change due to the spatial edit.
        profile = src.profile.copy()
        profile.update(
            width     = out_width,
            height    = out_height,
            transform = out_transform,
            tiled     = True,           # tile layout for efficient random access
            compress  = "deflate",      # lossless; swap for "lzw" if preferred
            bigtiff   = "IF_SAFER",     # auto BigTIFF when output exceeds 4 GB
        )

        total_chunks = math.ceil(src_height / chunk_rows)

        # ------------------------------------------------------------------
        # Write pixel data, tags, and colormap
        # ------------------------------------------------------------------
        with rasterio.open(output_path, "w", **profile) as dst:

            # Tags must be written while the dataset is open.
            if dataset_tags:
                dst.update_tags(**dataset_tags)
            if band_tags:
                dst.update_tags(1, **band_tags)

            for chunk_idx in range(total_chunks):
                src_row_start = chunk_idx * chunk_rows
                src_row_end   = min(src_row_start + chunk_rows, src_height)
                src_row_count = src_row_end - src_row_start

                # Read the full-width strip from the source.
                # Array shape: (1, src_row_count, src_width)
                src_data = src.read(
                    window=rasterio.windows.Window(
                        col_off=0,
                        row_off=src_row_start,
                        width=src_width,
                        height=src_row_count,
                    )
                )

                # Determine which rows in this strip survive and their
                # destination positions in the output file.
                for src_buf_offset, out_row, n_rows in out_row_segments(
                    src_row_start, src_row_count, src_height
                ):
                    # Extract the surviving row sub-strip from the chunk buffer.
                    # Shape: (1, n_rows, src_width)
                    sub = src_data[:, src_buf_offset: src_buf_offset + n_rows, :]

                    # Zero-fill the output buffer (unmapped columns stay 0 / nodata).
                    out_buf = np.zeros((1, n_rows, out_width), dtype=dtype)

                    # Fill Region A (left block) and Region B (shifted right block)
                    # into the output buffer in a single pass.
                    for src_col, count, out_off in src_cols_for_out_block(
                        0, out_width, src_width
                    ):
                        out_buf[:, :, out_off: out_off + count] = \
                            sub[:, :, src_col: src_col + count]

                    # Write the assembled strip to the correct output window.
                    dst.write(
                        out_buf,
                        window=rasterio.windows.Window(
                            col_off=0,
                            row_off=out_row,
                            width=out_width,
                            height=n_rows,
                        ),
                    )

                pct = (chunk_idx + 1) / total_chunks * 100
                print(
                    f"\r  Writing pixels : {pct:5.1f}%  "
                    f"(chunk {chunk_idx + 1}/{total_chunks})",
                    end="",
                    flush=True,
                )

            # ColorMap must be written before the dataset is closed.
            # rasterio silently discards write_colormap calls made after close().
            if colormap:
                dst.write_colormap(1, colormap)
                print("\n  ColorMap written.")

    # RAT cannot be written through rasterio's public API.
    # We reopen the now-closed output file with GDAL in update mode and
    # attach the RAT there.  This is safe because rasterio has already
    # released all file handles.
    print("  Copying RAT ...")
    copy_rat(input_path, output_path, band=1)

    print(f"Done.  Output file: {output_path}")


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Crop and shift a single-band GeoTIFF by removing a contiguous "
            "column gap and a contiguous row gap, then closing the gaps by "
            "shifting the right/bottom portions inward."
        )
    )
    parser.add_argument("input",  help="Path to the source GeoTIFF file.")
    parser.add_argument("output", help="Path for the processed output GeoTIFF.")
    parser.add_argument(
        "--chunk-rows",
        type=int,
        default=512,
        metavar="N",
        help=(
            "Number of source rows to process per iteration (default: 512). "
            "Decrease this value if RAM is limited; increase it for faster I/O "
            "on systems with ample memory."
        ),
    )
    args = parser.parse_args()
    process(args.input, args.output, chunk_rows=args.chunk_rows)


if __name__ == "__main__":
    main()
