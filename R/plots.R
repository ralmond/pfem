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
    results <- dplyr::filter(results,method==input$method,
                      subj==input$subj)
    summary <- dplyr::filter(summary,method%in%input$traces,
                      subj==input$subj)
    ggplot2::ggplot(results, ggplot2::aes(x=time,y=theta,colour=weights)) +
      ggplot2::geom_point() +
      ggplot2::scale_colour_gradient(low="gray90",high="gray10") +
      ggplot2::geom_line(data=summary,
                mapping=aes(x=time,y=theta_bar,linetype=method),
                color="blue")+
      ggplot2::geom_line(data=summary,
                mapping=aes(x=time,y=theta0),
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
      ggplot2::ggplot(biastab,ggplot2::aes(y=scaledBias,x=method,colour=subj)) +
        ggplot2::scale_y_continuous(limits=c(-3,3)) +
        ggplot2::geom_point() +
        ggplot2::geom_line(data=dplyr::filter(biastab,subj%in%input$subjcc),
                  mapping=ggplot2::aes(x=as.numeric(method),
                              y=scaledBias,colour=subj)) +
        ggplot2::geom_hline(ggplot2::aes(yintercept=1))
    )

    output$table <- shiny::renderTable(
      biastab |> dplyr::group_by(method) |>
      dplyr::summarize(bias=mean(bias),mse=mean(bias^2),
                t=mean(scaledBias), t2 = sum(scaledBias^2)) |>
      round(3))
  }
  shiny::shinyApp(ui,server)
}
