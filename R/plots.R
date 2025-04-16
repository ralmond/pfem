bySubjPlot <- function(results,theta0=NULL,
                       lowcol="gray90",highcol="gray10",
                       aveshape=4,avecol="red",avesize=2,
                       trueshape=8,truecol="blue",truesize=2) {
  summar <- avePart(results)
  if (!missing(theta0)) summar$theta0 <- theta0
  plot <-
    ggplot2::ggplot(results,
                    ggplot2::aes(x=subj,y=theta,colour=weights,
                                 frame=occ)) +
    ggplot2::geom_point() +
    ggplot2::scale_colour_gradient(low=lowcol,high=highcol) +
    ggplot2::geom_point(data=summar,
                        mapping=ggplot2::aes(x=subj,y=theta_bar,frame=occ),
                        shape=aveshape,col=avecol,size=avesize)
  if (!missing(theta0))
    plot <- plot +
      ggplot2::geom_point(data=summar,
                          mapping=ggplot2::aes(x=subj,y=theta0,frame=occ),
                          shape=trueshape,col=truecol,size=truesize)
  plot
}

byTimePlot <- function(results,theta0=NULL,
                       lowcol="gray90",highcol="gray10",
                       avetype=2,avecol="red",linewidth=1,
                       truetype=1,truecol="blue") {
  summar <- avePart(results)
  if (!missing(theta0)) summar$theta0 <- theta0
  plot <-
    ggplot2::ggplot(results,
                    ggplot2::aes(x=time,y=theta,colour=weights,
                                 frame=subj)) +
    ggplot2::geom_point() +
    ggplot2::scale_colour_gradient(low=lowcol,high=highcol) +
    ggplot2::geom_line(data=summar,
                        mapping=ggplot2::aes(x=time,y=theta_bar,frame=subj),
                        linetype=avetype,col=avecol,linewidth=linewidth)
  if (!missing(theta0))
    plot <- plot +
      ggplot2::geom_line(data=summar,
                          mapping=ggplot2::aes(x=time,y=theta0,frame=subj),
                          linetype=truetype,col=truecol,linewidth=linewidth)
  plot
}
