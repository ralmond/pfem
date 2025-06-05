bySubjPlot <- function(results,theta0=NULL,
                       lowcol="gray90",highcol="gray10",
                       aveshape=4,avecol="red",avesize=2,
                       trueshape=8,truecol="blue",truesize=2) {
  summar <- avePart(results)
  if (!missing(theta0)) summar$theta0 <- theta0
  plot <-
    ggplot2::ggplot(results,
                    ggplot2::aes(x=results$subj,y=results$theta,
                                 colour=results$weights,
                                 frame=results$occ)) +
    ggplot2::geom_point() +
    ggplot2::scale_colour_gradient(low=lowcol,high=highcol) +
    ggplot2::geom_point(data=summar,
                        mapping=ggplot2::aes(x=summar$subj,
                                             y=summar$theta_bar,
                                             frame=summar$occ),
                        shape=aveshape,col=avecol,size=avesize)
  if (!missing(theta0))
    plot <- plot +
      ggplot2::geom_point(data=summar,
                          mapping=ggplot2::aes(x=summar$subj,
                                               y=summar$theta0,
                                               frame=summar$occ),
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
                    ggplot2::aes(x=results$time,y=results$theta,
                                 colour=results$weights,
                                 frame=results$subj)) +
    ggplot2::geom_point() +
    ggplot2::scale_colour_gradient(low=lowcol,high=highcol) +
    ggplot2::geom_line(data=summar,
                        mapping=ggplot2::aes(x=summar$time,
                                             y=summar$theta_bar,
                                             frame=summar$subj),
                        linetype=avetype,col=avecol,linewidth=linewidth)
  if (!missing(theta0))
    plot <- plot +
      ggplot2::geom_line(data=summar,
                          mapping=ggplot2::aes(x=summar$time,y=summar$theta0,
                                               frame=summar$subj),
                          linetype=truetype,col=truecol,linewidth=linewidth)
  plot
}

methodTimeGadget <- function(results,summary,methods) {

  subjlist <- unique(results$subj)
  ui <- shiny::fillPage(
    shiny::inputPanel(
      shiny::selectInput("method","Method: ",methods,
                  selected=methods[1], multiple=FALSE),
      shiny::selectInput("subj","Subject:",subjlist,subjlist[1]),
      shiny::selectInput("traces","Selected Traces: ",methods,
                  selected=methods[1], multiple=TRUE,selectize=FALSE)
    ),
    shiny::plotOutput("plot")
  )
  server <- function(input,output,session) {
    output$plot <-  shiny::renderPlot({
    results <- dplyr::filter(results,results$method==input$method,
                      results$subj==input$subj)
    summary <- dplyr::filter(summary,summary$method%in%input$traces,
                      summary$subj==input$subj)
    ggplot2::ggplot(results,
                    ggplot2::aes(x=results$time,y=results$theta,
                                 colour=results$weights)) +
      ggplot2::geom_point() +
      ggplot2::scale_colour_gradient(low="gray90",high="gray10") +
      ggplot2::geom_line(data=summary,
                mapping=ggplot2::aes(x=summary$time,y=summary$theta_bar,
                                     linetype=summary$method),
                color="blue")+
      ggplot2::geom_line(data=summary,
                mapping=ggplot2::aes(x=summary$time,y=summary$theta0),
                color="red")+
     ggplot2::guides(colour="none")
    })
  }
  shiny::shinyApp(ui,server)
}


biasDisplayGadget <- function(biastab) {
  subjlist <- unique(biastab$subj)

  ui <- shiny::fillPage(
    shiny::inputPanel(
      shiny::selectInput("subj","Subject:",1:10,multiple=TRUE,selectize=FALSE)
    ),
    shiny::plotOutput("plot"),shiny::tableOutput("table")
  )

  server <- function(input,output,session) {
    output$plot <- shiny::renderPlot(
      ggplot2::ggplot(biastab,ggplot2::aes(y=biastab$scaledBias,
                                           x=biastab$method,
                                           colour=biastab$subj)) +
        ggplot2::scale_y_continuous(limits=c(-3,3)) +
        ggplot2::geom_point() +
        ggplot2::geom_line(data=dplyr::filter(biastab,
                                              biastab$subj%in%input$subjcc),
                  mapping=ggplot2::aes(x=as.numeric(.data$method),
                              y=.data$scaledBias,colour=.data$subj)) +
        ggplot2::geom_hline(ggplot2::aes(yintercept=0))
    )

    output$table <- shiny::renderTable(
      biastab |> dplyr::group_by(.data$method) |>
      dplyr::summarize(bias=mean(.data$bias),mse=mean(.data$bias^2),
                t=mean(.data$scaledBias), t2 = sum(.data$scaledBias^2)) |>
      round(3))
  }
  shiny::shinyApp(ui,server)
}
