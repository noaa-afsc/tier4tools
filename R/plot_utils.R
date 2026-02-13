#' Tier 4 default ggplot theme
#'
#' A SAFE-friendly theme with minimal gridlines and consistent
#' strips/legend.
#'
#' @return A ggplot2 theme object.
#' @export
theme_tier4 <- function(base_size = 12) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      legend.position = "right",
      # Increase legend key size
      legend.key.size = grid::unit(1.2, "lines"),
      legend.spacing.y = grid::unit(0.4, "lines")
    )
}

#' Tier 4 default colorblind-friendly palette (ggthemes "colorblind" hex codes)
#'
#' Discrete palette intended for species, scenarios, or other categorical
#' aesthetics.
#'
#' @param n Number of colors requested.
#' @param reverse Logical. If TRUE, reverse the palette order.
#' @param recycle Logical. If TRUE and n exceeds palette length, recycle colors
#'   (default TRUE).
#' @param warn Logical. If TRUE, warn when recycling is needed (default TRUE).
#'
#' @return Character vector of hex color codes.
#' @keywords internal
tier4_palette_cb <- function(
    n,
    reverse = FALSE,
    recycle = TRUE,
    warn = TRUE
) {
  pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  if (reverse) pal <- rev(pal)

  if (!is.numeric(n) || length(n) != 1 || !is.finite(n) || n < 1) {
    stop("n must be a single positive number.")
  }
  n <- as.integer(n)

  if (n <= length(pal)) return(pal[seq_len(n)])

  if (!isTRUE(recycle)) {
    stop("Requested n = ", n, " colors, but tier4 palette has only ", length(pal), ".")
  }
  if (isTRUE(warn)) {
    warning("tier4_palette_cb(): n = ", n, " exceeds palette length (", length(pal),
            "); recycling colors.")
  }
  rep(pal, length.out = n)
}

#' Tier 4 discrete color scale
#'
#' Uses the tier4 colorblind-friendly discrete palette (ggthemes colorblind hex codes).
#'
#' @param ...,reverse,recycle,warn Passed to internal palette helper.
#' @param na.value Color for NA values.
#'
#' @return A ggplot2 scale.
#' @export
scale_color_tier4 <- function(..., reverse = FALSE, recycle = TRUE, warn = TRUE, na.value = "grey50") {
  ggplot2::discrete_scale(
    aesthetics = "colour",
    scale_name = "tier4_cb",
    palette = function(n) tier4_palette_cb(n, reverse = reverse, recycle = recycle, warn = warn),
    na.value = na.value,
    ...
  )
}

#' Tier 4 discrete fill scale
#'
#' Uses the tier4 colorblind-friendly discrete palette (ggthemes colorblind hex
#' codes).
#'
#' @param ...,reverse,recycle,warn Passed to internal palette helper.
#' @param na.value Fill color for NA values.
#'
#' @return A ggplot2 scale.
#' @export
scale_fill_tier4 <- function(..., reverse = FALSE, recycle = TRUE, warn = TRUE, na.value = "grey50") {
  ggplot2::discrete_scale(
    aesthetics = "fill",
    scale_name = "tier4_cb",
    palette = function(n) tier4_palette_cb(n, reverse = reverse, recycle = recycle, warn = warn),
    na.value = na.value,
    ...
  )
}
