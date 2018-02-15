# Functions below are used to visualize the output of the ALStructure.

#' Generate a scatterplot of rows of \eqn{\boldsymbol{Q}}{Q} and \eqn{\boldsymbol{\hat{Q}}}{Q_hat}.
#'
#' Generate a scatterplot of rows of \eqn{\boldsymbol{Q}}{Q} and \eqn{\boldsymbol{\hat{Q}}}{Q_hat}
#' for a visual comparison.
#' This can be used to compare the structure of \eqn{\boldsymbol{\hat{Q}}}{Q_hat} fitted using
#' \code{ALStructure} to the truth (in the case of simulated data), or it can be
#' used to compare the outputs of two fitting algorithms. The plots are written
#' into the file specified by \code{output_file}.
#'
#' @param Q The "true" \eqn{d \times n}{d x n} underlying admixture matrix
#' @param Q_hat The estimated \eqn{k \times n}{d x n} admixture matrix to compare to Q
#' @param output_file The file to which the generated plot will be written
#' @param dims A vector of length 2 that specifies which 2 rows will be plotted
#' @param res The resolution of the output figure
#' @param rand_held_out The number of randomly selected points which will be
#'        colored differently from the others. This feather is to allow for a
#'        clear visual demonstration of how well individual points concur
#'        between \eqn{\boldsymbol{Q}}{Q} and \eqn{\boldsymbol{\hat{Q}}}{Q_hat}.
#'
#' @return A \code{ggplot} object \code{p_all} containing the scatterplot of both the estimate
#'         and ground truth side by side.
#'
#' @keywords internal
compare_Q <- function(Q, Q_hat, output_file = NaN, dims = c(1,2), res = 100,
                      rand_held_out = 0){
  library(gridExtra)
  ggplot2::theme_set(ggplot2::theme_bw())
  axis_x <- sprintf("Q%d", dims[1]);
  axis_y <- sprintf("Q%d", dims[2]);
  aes_x <- sprintf("X%d", dims[1]);
  aes_y <- sprintf("X%d", dims[2]);

  Q <- data.frame(t(Q)); Q_hat <- data.frame(t(Q_hat))
  n <- dim(Q)[1]

  if (rand_held_out == 0){
  p1 <- ggplot2::ggplot(data = Q) +
        ggplot2::geom_point(mapping = ggplot2::aes_string(x = aes_x, y = aes_y),
                            color="tomato", alpha=0.3) +
        ggplot2::xlim(0,1) + ggplot2::ylim(0,1) +
        ggplot2::xlab(axis_x) + ggplot2::ylab(axis_y) +
        ggplot2::ggtitle("Truth") + ggplot2::coord_fixed(ratio = 1)

  p2 <- ggplot2::ggplot(data = Q_hat) +
        ggplot2::geom_point(mapping = ggplot2::aes_string(x = aes_x, y = aes_y),
                            color="tomato", alpha=0.3) +
        ggplot2::xlim(0,1) + ggplot2::ylim(0,1) +
        ggplot2::xlab(axis_x) + ggplot2::ylab(axis_y) +
        ggplot2::ggtitle("Fit") + ggplot2::coord_fixed(ratio = 1)

  } else {
    s <- sample(n, rand_held_out)
    held_out <- vector(mode = "logical", length = n); held_out[s] <- TRUE
    held_hout <- as.factor(held_out)

    p1 <- ggplot2::ggplot(data = Q) +
          ggplot2::geom_point(mapping = ggplot2::aes_string(x = aes_x, y = aes_y,
                              color = "held_out", alpha = "held_out")) +
          ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
          ggplot2::xlab(axis_x) + ggplot2::ylab(axis_y) +
          ggplot2::ggtitle("Truth") + ggplot2::theme(legend.position="none") +
          ggplot2::coord_fixed(ratio = 1)

    p1 <- p1 + ggplot2::scale_color_manual(breaks = c("TRUE", "FALSE"), values = c("tomato", "slateblue")) +
          ggplot2::scale_alpha_manual(breaks = c("TRUE", "FALSE"), values = c(0.3, 1))



    p2 <- ggplot2::ggplot(data = Q_hat) +
          ggplot2::geom_point(mapping = ggplot2::aes_string(x = aes_x, y = aes_y,
                              color = "held_out", alpha = "held_out")) +
          ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
          ggplot2::xlab(axis_x) + ggplot2::ylab(axis_y) +
          ggplot2::ggtitle("Fit") + ggplot2::theme(legend.position="none") +
          ggplot2::coord_fixed(ratio = 1)

    p2 <- p2 + ggplot2::scale_color_manual(breaks = c("TRUE", "FALSE"), values = c("tomato", "slateblue")) +
          ggplot2::scale_alpha_manual(breaks = c("TRUE", "FALSE"), values = c(0.3, 1))
  }

  p_all <- grid.arrange(p1, p2, ncol = 2)

  if(!is.nan(output_file)){
    png(output_file, width = 4.2, height = 4, units = 'in', res = res)
    grid.arrange(p1, p2, ncol = 2)
    dev.off()
  }

  p_all
}

#' Multiple plot function
#'
#' Plot a colection of ggplot2 plot objects into a simple grid layout with a
#' specified number of columns
#'
#' @param ... the ggplot2 plot objects to be plotted
#' @param cols the number of columns into which the plots will be arranged
#'
#' @return a ggplot2 object with supplied plots arranged in a grid
#'
#' @keywords internal

multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
      ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
        layout.pos.col = matchidx$col))
    }
  }
}
