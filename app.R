#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(BiocManager)
options(repos = BiocManager::repositories())

library(shiny)
library(igraph)
library(tidygraph)
library(ggraph)
library(scran)

load("D_exp_z.RData")
load("D_spa_z.RData")
load("X1_z.RData")
load("X2_z.RData")
load("top10_var_counts.RData")
load("centroids.RData")
var_genes <- colnames(var_counts)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Mouse Cortex scRNA-seq + Spatial Vizualization"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("w",
                        "Mixing weight",
                        min = 0,
                        max = 1,
                        value = 0),
            radioButtons("feature",
                         "Display Feature",
                         choices = c("Cell Clusters",
                                     "Variable Genes")),
            selectInput("gene",
                        "Variable Genes",
                        choices = var_genes,
                        selected = "Mag"),
            hr(),
            numericInput("knn",
                         "Number of nearest neighbors",
                         min = 3,
                         max = floor(sqrt(nrow(D1))),
                         value = floor((sqrt(nrow(D1)) + 3)/2))
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot",
                      brush = brushOpts(id = "brush_pts")),
           plotOutput("spatPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
        
        w = input$w
        D = (1-w)*D1 + w*D2
        G <- buildKNNGraph(D, k = input$knn)
        XY <- (1-w)*X1_z + (w)*X2_z
        
        ingene = input$gene
        incounts = var_counts[,ingene]
        
        fit_l <- cluster_louvain(G)
        clusts_l <- as.factor(fit_l$membership)
        g <- G %>%
            as_tbl_graph() %>%
            activate(nodes) %>%
            mutate(count = incounts,
                   clust = clusts_l)
        XY_df <- as.data.frame(XY)
        XY_df$clust <- clusts_l
        XY_df$count <- incounts
        colnames(XY_df) <- c("X","Y","clust","count")
        if(input$feature == "Cell Clusters")
        {
            p <- ggraph(g,layout = "kk") + 
                geom_node_point(aes(color = clust)) +
                ggtitle("Blended cell distances") + 
                theme(legend.position = "none")   
        }
        if(input$feature == "Variable Genes")
        {
            p <- ggraph(g,layout = "kk") + 
                geom_node_point(aes(color = count)) +
                ggtitle("Blended cell distances") + 
                theme(legend.position = "none") 
        }
        
        pg <- ggplot_build(p)
        coords <- pg[[1]][[1]][,c("x","y")]
        save(coords,file = "coords.RData")
        p
    })
    
    output$spatPlot <- renderPlot({
        xs <- c(input$brush_pts$xmin,input$brush_pts$xmax)
        ys <- c(input$brush_pts$ymin,input$brush_pts$ymax)
        load("coords.RData")
        inds_sel <- coords$x > xs[1] & coords$x < xs[2] & coords$y > ys[1] & coords$y < ys[2]
        cents_sel <- centroids[inds_sel,]
        ggplot() + 
            geom_point(data = centroids, aes(x = X, y = Y), alpha = 0.1) +
            geom_point(data = cents_sel, aes(x = X, y = Y)) +
            ggtitle("Selected cell locations") +
            xlab(NULL) + ylab(NULL) +
            theme(axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank())
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
