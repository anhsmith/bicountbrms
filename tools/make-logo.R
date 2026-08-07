# Generate the bicountbrms hex sticker: man/figures/logo.svg and logo.png.
# Run from the package root.
#
# The motif is the model itself. Counts come in pairs. Within a pair the grey
# blocks are the shared latent component -- equal in both columns by
# construction. Above them sit the two private components, blue for the first
# source and orange for the second; their difference is the disagreement the
# families estimate.

W <- 519.6   # hex width  (1.732 : 2, the hexb.in ratio)
H <- 600     # hex height

col_shared <- "#B4AA98"  # warm grey -- shared latent component
col_src1   <- "#5BA5C7"  # blue      -- source 1, private
col_src2   <- "#E28A45"  # orange    -- source 2, private
col_ink    <- "#2F5F7A"
col_bg     <- "#FCFAF6"

hex_pts <- paste(
  sprintf(
    "%.2f,%.2f",
    c(W / 2, W, W, W / 2, 0, 0),
    c(0, H / 4, 3 * H / 4, H, 3 * H / 4, H / 4)
  ),
  collapse = " "
)

# ---- the pairs --------------------------------------------------------------
# one row per pair: shared blocks, then the private blocks of each source

pairs <- list(
  c(shared = 1, src1 = 1, src2 = 3),
  c(shared = 4, src1 = 2, src2 = 1),
  c(shared = 2, src1 = 3, src2 = 2)
)

sq        <- 40                    # block side
gap_y     <- 9
gap_in    <- 8                     # gap within a pair
gap_out   <- 40                    # gap between pairs
pitch_y   <- sq + gap_y
pair_w    <- 2 * sq + gap_in
pitch_x   <- pair_w + gap_out

base_y <- 385                                        # common baseline
span   <- length(pairs) * pair_w + (length(pairs) - 1) * gap_out
x0     <- (W - span) / 2

blocks <- character(0)

blk <- function(x, i, fill) {
  sprintf(
    '    <rect x="%.2f" y="%.2f" width="%d" height="%d" rx="7" fill="%s"/>',
    x, base_y - i * pitch_y + gap_y, sq, sq, fill
  )
}

for (p in seq_along(pairs)) {
  s  <- pairs[[p]]
  xp <- x0 + (p - 1) * pitch_x

  for (m in 0:1) {
    x <- xp + m * (sq + gap_in)
    for (i in seq_len(s[["shared"]])) {
      blocks <- c(blocks, blk(x, i, col_shared))
    }
    private <- if (m == 0) s[["src1"]] else s[["src2"]]
    for (i in seq_len(private)) {
      blocks <- c(blocks, blk(x, s[["shared"]] + i, if (m == 0) col_src1 else col_src2))
    }
  }
}

# ---- assemble ---------------------------------------------------------------

svg <- c(
  sprintf(
    '<svg xmlns="http://www.w3.org/2000/svg" width="%.0f" height="%d" viewBox="0 0 %.2f %d">',
    W, H, W, H
  ),
  '  <defs>',
  sprintf('    <clipPath id="hex"><polygon points="%s"/></clipPath>', hex_pts),
  '  </defs>',
  sprintf('  <polygon points="%s" fill="%s"/>', hex_pts, col_bg),
  '  <g clip-path="url(#hex)">',
  blocks,
  sprintf(
    '    <text x="%.2f" y="450" text-anchor="middle" font-family="Helvetica, Arial, sans-serif" font-size="58" font-weight="600" letter-spacing="0.5" fill="%s">bicountbrms</text>',
    W / 2, col_ink
  ),
  '  </g>',
  sprintf(
    '  <polygon points="%s" fill="none" stroke="%s" stroke-width="20" clip-path="url(#hex)"/>',
    hex_pts, col_ink
  ),
  '</svg>'
)

dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)
writeLines(svg, "man/figures/logo.svg")
rsvg::rsvg_png("man/figures/logo.svg", "man/figures/logo.png", width = 520, height = 600)
cat("written\n")
