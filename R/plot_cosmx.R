#' Plot spatial expression of multiple genes in CosMx data
#'
#' Creates spatial transcriptomics visualizations for CosMx data, displaying
#' the expression of multiple genes as colored points overlaid on cell polygons.
#' The function generates both a full view of the entire dataset and an optional
#' zoomed-in view of a specified region.
#'
#' @param polygons A data frame containing cell polygon coordinates. Must include
#'   columns: \code{x_global_px}, \code{y_global_px}, \code{fov}, and \code{cellID}.
#' @param tx A data frame containing transcript location data. Must include
#'   columns: \code{x_global_px}, \code{y_global_px}, and \code{target} (gene name).
#' @param genes A character vector of gene names to plot. Each gene will be
#'   displayed in a different color.
#' @param colors Optional character vector of colors for each gene. If \code{NULL},
#'   colors are automatically generated using \code{rainbow()}. Must match the
#'   length of \code{genes}.
#' @param zoom.x Optional numeric value in [0, 1] specifying the relative x position
#'   of the zoom window center (0 = left, 1 = right, 0.5 = center). Default: 0.5.
#' @param zoom.y Optional numeric value in [0, 1] specifying the relative y position
#'   of the zoom window center (0 = top, 1 = bottom, 0.5 = center). Default: 0.5.
#' @param zoom.width.x Numeric. Width of the zoom window in pixels (x-direction).
#'   Default: 1000.
#' @param zoom.width.y Numeric. Height of the zoom window in pixels (y-direction).
#'   Default: 1000.
#' @param microns_per_pixel Numeric. Conversion factor from pixels to microns.
#'   Default: 0.12028.
#' @param scale_length_um_full Numeric. Length of the scale bar in microns for
#'   the full view. Default: 50.
#' @param scale_length_um_zoom Numeric. Length of the scale bar in microns for
#'   the zoom view. Default: 25.
#' @param label_size_full Numeric. Font size for the scale bar label in the full
#'   view. Default: 22.
#' @param label_size_zoom Numeric. Font size for the scale bar label in the zoom
#'   view. Default: 25.
#' @param point_size_full Numeric. Size of gene expression points in the full view.
#'   Default: 0.3.
#' @param point_size_zoom Numeric. Size of gene expression points in the zoom view.
#'   Default: 2.
#' @param polygon_linewidth_full Numeric. Line width for cell polygon borders in
#'   the full view. Default: 0.3.
#' @param polygon_linewidth_zoom Numeric. Line width for cell polygon borders in
#'   the zoom view. Default: 1.
#' @param show_zoom_rect Logical. Whether to display a rectangle on the full view
#'   indicating the zoom region. Default: \code{TRUE}.
#' @param zoom.rect.color Character. Color of the zoom rectangle border.
#'   Default: "white".
#' @param zoom.rect.alpha Numeric. Transparency of the zoom rectangle fill (0-1).
#'   Default: 0 (fully transparent).
#' @param zoom.rect.linewidth Numeric. Line width of the zoom rectangle border.
#'   Default: 1.
#' @param zoom.rect.linetype Character. Line type of the zoom rectangle border
#'   (e.g., "solid", "dashed", "dotted"). Default: "dashed".
#'
#' @return A list containing two ggplot2 objects:
#'   \itemize{
#'     \item \code{full}: A ggplot2 object showing the full spatial view with all
#'       genes plotted as colored points.
#'     \item \code{zoom}: A ggplot2 object showing the zoomed-in view, or \code{NULL}
#'       if zoom parameters are not provided.
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' plots <- plotCosMxGenes(
#'     polygons = cell_polygons,
#'     tx = transcript_data,
#'     genes = c("Gene1", "Gene2", "Gene3"),
#'     colors = c("red", "green", "blue"),
#'     zoom.x = 0.5,
#'     zoom.y = 0.5
#' )
#'
#' # Display the plots
#' plots$full
#' plots$zoom
#' }
plotCosMxGenes <- function(
    polygons,
    tx,
    genes,
    colors = NULL,
    zoom.x = NULL,
    zoom.y = NULL,
    zoom.width.x = 1000,
    zoom.width.y = 1000,
    microns_per_pixel = 0.12028,
    scale_length_um_full = 50,
    scale_length_um_zoom = 25,
    label_size_full = 22,
    label_size_zoom = 25,
    point_size_full = 0.3,
    point_size_zoom = 2,
    polygon_linewidth_full = 0.3,
    polygon_linewidth_zoom = 1,
    show_zoom_rect = TRUE,
    zoom.rect.color = "white",
    zoom.rect.alpha = 0,
    zoom.rect.linewidth = 1,
    zoom.rect.linetype = "dashed") {
    # Default colors if not provided
    if (is.null(colors)) {
        colors <- rainbow(length(genes))
    }

    # Check that colors match genes length
    if (length(colors) != length(genes)) {
        stop("Length of colors must match length of genes")
    }

    # Filter tx data for each gene
    tx.data_list <- lapply(genes, function(gene) {
        subset(tx, target == gene)
    })
    names(tx.data_list) <- genes

    # Calculate full view range
    x_range <- range(polygons$x_global_px, na.rm = TRUE)
    y_range <- range(polygons$y_global_px, na.rm = TRUE)
    x_min <- x_range[1]
    x_max <- x_range[2]
    y_min <- y_range[1]
    y_max <- y_range[2]
    x_span <- x_max - x_min
    y_span <- y_max - y_min

    # Set default zoom coordinates to center (0.5, 0.5) if not provided
    if (is.null(zoom.x)) {
        zoom.x <- 0.5
    }
    if (is.null(zoom.y)) {
        zoom.y <- 0.5
    }

    # Validate zoom coordinates are in [0, 1]
    if (zoom.x < 0 || zoom.x > 1 || zoom.y < 0 || zoom.y > 1) {
        stop("zoom.x and zoom.y must be in the range [0, 1]")
    }

    # Convert relative position [0,1] to absolute pixel coordinates
    # zoom.x/y represents the relative position of the zoom window center
    zoom.x.abs <- x_min + zoom.x * x_span
    zoom.y.abs <- y_min + zoom.y * y_span

    # Calculate zoom window with center at zoom.x/y, ensuring it stays within bounds
    zoom.x.start <- zoom.x.abs - zoom.width.x / 2
    zoom.y.start <- zoom.y.abs - zoom.width.y / 2

    # Clamp zoom window to stay within full view bounds
    if (zoom.x.start < x_min) {
        zoom.x.start <- x_min
    }
    if (zoom.x.start + zoom.width.x > x_max) {
        zoom.x.start <- x_max - zoom.width.x
    }
    if (zoom.y.start < y_min) {
        zoom.y.start <- y_min
    }
    if (zoom.y.start + zoom.width.y > y_max) {
        zoom.y.start <- y_max - zoom.width.y
    }

    # Create zoom window
    zoom.window <- c(zoom.x.start, zoom.x.start + zoom.width.x, zoom.y.start, zoom.y.start + zoom.width.y)

    # Create scale bar for full view
    scale.bar.full <- sputilsR::make_scale_bar(
        x_vals = polygons$x_global_px,
        y_vals = polygons$y_global_px,
        microns_per_pixel = microns_per_pixel,
        scale_length_um = scale_length_um_full,
        label_size = label_size_full
    )

    # Build full view plot
    p_full <- ggplot2::ggplot() +
        ggplot2::geom_polygon(
            data = polygons,
            mapping = ggplot2::aes(x = x_global_px, y = y_global_px, group = interaction(fov, cellID)),
            fill = "white",
            color = "darkgrey",
            linewidth = polygon_linewidth_full,
            alpha = 0
        )

    # Add points for each gene
    for (i in seq_along(genes)) {
        p_full <- p_full +
            ggplot2::geom_point(
                data = tx.data_list[[i]],
                mapping = ggplot2::aes(x = x_global_px, y = y_global_px),
                color = colors[i],
                size = point_size_full
            )
    }

    # Add zoom rectangle if specified
    if (show_zoom_rect && !is.null(zoom.window)) {
        p_full <- p_full +
            ggplot2::annotate(
                "rect",
                xmin = zoom.window[1],
                xmax = zoom.window[2],
                ymin = zoom.window[3],
                ymax = zoom.window[4],
                color = zoom.rect.color,
                alpha = zoom.rect.alpha,
                linewidth = zoom.rect.linewidth,
                linetype = zoom.rect.linetype
            )
    }

    p_full <- p_full +
        scale.bar.full$bg + scale.bar.full$rect + scale.bar.full$label +
        ggplot2::scale_fill_viridis_c() +
        ggplot2::coord_equal() +
        ggplot2::theme_void() +
        ggplot2::theme(panel.background = ggplot2::element_rect(fill = "black")) +
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggplot2::scale_y_continuous(expand = c(0, 0))

    # Create zoom view if zoom window is provided
    p_zoom <- NULL
    if (!is.null(zoom.window)) {
        # Filter polygons and tx data for zoom window
        zoom.polygons <- subset(
            polygons,
            x_global_px > zoom.window[1] & x_global_px < zoom.window[2] &
                y_global_px > zoom.window[3] & y_global_px < zoom.window[4]
        )
        zoom.polygons <- zoom.polygons %>%
            dplyr::group_by(cellID) %>%
            dplyr::filter(dplyr::n() > 5)

        zoom.tx.data_list <- lapply(tx.data_list, function(tx.data) {
            subset(
                tx.data,
                x_global_px > zoom.window[1] & x_global_px < zoom.window[2] &
                    y_global_px > zoom.window[3] & y_global_px < zoom.window[4]
            )
        })

        # Create scale bar for zoom view
        scale.bar.zoom <- sputilsR::make_scale_bar(
            x_vals = zoom.polygons$x_global_px,
            y_vals = zoom.polygons$y_global_px,
            microns_per_pixel = microns_per_pixel,
            scale_length_um = scale_length_um_zoom,
            label_size = label_size_zoom
        )

        # Build zoom view plot
        p_zoom <- ggplot2::ggplot() +
            ggplot2::geom_polygon(
                data = zoom.polygons,
                mapping = ggplot2::aes(x = x_global_px, y = y_global_px, group = interaction(fov, cellID)),
                fill = "white",
                color = "darkgrey",
                linewidth = polygon_linewidth_zoom,
                alpha = 0
            )

        # Add points for each gene
        for (i in seq_along(genes)) {
            p_zoom <- p_zoom +
                ggplot2::geom_point(
                    data = zoom.tx.data_list[[i]],
                    mapping = ggplot2::aes(x = x_global_px, y = y_global_px),
                    color = colors[i],
                    size = point_size_zoom
                )
        }

        p_zoom <- p_zoom +
            scale.bar.zoom$bg + scale.bar.zoom$rect + scale.bar.zoom$label +
            ggplot2::scale_fill_viridis_c() +
            ggplot2::coord_equal() +
            ggplot2::theme_void() +
            ggplot2::theme(panel.background = ggplot2::element_rect(fill = "black")) +
            ggplot2::scale_x_continuous(expand = c(0, 0)) +
            ggplot2::scale_y_continuous(expand = c(0, 0))
    }

    # Return list of plots
    return(list(full = p_full, zoom = p_zoom))
}
