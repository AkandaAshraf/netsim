#' This fucntion allowes you to generate three different types of networks with same number of vertices and edges
#' @param ScaleFreePowerRange A vector containing powers for scale free networks, e.g. c(1.5,1.75)
#' @param SmallWorldProbability A vector containing probabilities for small world networks, e.g. c(0.5,0.6)
#' @param VerticesVector A vector containing number of vertices to be used to generate all three network types(please note that number of edges will be calculated to match number of edges and verices for all three network types, see documentation for further information), e.g. c(10,100)
#' @param SampleSize Integer number defining how many sample will be generated for each network types
#' @param edgesEachScaleFreeInteger number defining how many edges to be added in each step of scale free networks
#' @param savingDir A string with the directory where to save graph objects and plots, the directory folder need to be created beforehand
#' @param plotGraph if "y" then plots will also be created, only line plots to generate point plots use
#' @import igraph
#' @import ggplot2
#' @export

gen_graphs<- function(ScaleFreePowerRange,SmallWorldProbability,VerticesVector,SampleSize,edgesEachScaleFree,savingDir,plotGraph="y")
{
  #require('igraph')

  # Start the clock!
  ptm <- proc.time()

  print("generating Scalefree networks(preferential attachment, 2 edges each step) and calculating different properties")

  # EdgesVector = (VerticesVector * 2)-3

  EdgesVector = ((VerticesVector*edgesEachScaleFree)-((edgesEachScaleFree*(edgesEachScaleFree+1))/2))

  NeiSmallWorld = edgesEachScaleFree

  numberOFDelforSmallWorld = (VerticesVector*(edgesEachScaleFree-NeiSmallWorld))+((edgesEachScaleFree*(edgesEachScaleFree+1))/2)
  print(numberOFDelforSmallWorld)

  ##numberOFDelforSmallWorld = ((edgesEachScaleFree*(edgesEachScaleFree+1))/2)
  ##print(numberOFDelforSmallWorld)
  ####Generate Networks

  for(power in ScaleFreePowerRange)
  {
    ScaleFreeGraphsEdges = lapply(1:SampleSize,function(i) BAGraphScalefreePropertiesEdges(VectorVertices = VerticesVector,p1=power,p2=power,PstepSize=0.25,EdgesEachStep = edgesEachScaleFree))
    print(paste("Serializing and saving Network objects Object:Scalefree",power,sep = ""))
    saveRDS(ScaleFreeGraphsEdges, file = paste(savingDir,"Scalefree",power,sep = ""), ascii = FALSE, version = NULL,
            compress = TRUE, refhook = NULL)

  }

  print("generating and calculating different properties of random graph networks(gnm) with same number of edges as in Scalefree ")

  RandomGraphsEdges = lapply(1:SampleSize,function(i) ERGraphRandomPropertiesGNM2(VectorVertices = VerticesVector,VectorEdges = EdgesVector,PstepSize=1))
  print(paste("Serializing and saving Network objects Object: RandomGraphsEdges"))
  saveRDS(RandomGraphsEdges, file =paste(savingDir,"randomGraphEdgeComparisn",sep=""), ascii = FALSE, version = NULL,
          compress = TRUE, refhook = NULL)



  print("generating small World Networks and calculating different properties lattice Dimention-1, Nei-2, also deleting randomly selected edges to match with other two")

  for(Probability in SmallWorldProbability)
  {
    SmallWorldEdges = lapply(1:SampleSize,function(i) SmallWorldGraphPropertiesEdges(SizeVector = VerticesVector,latticeNei = NeiSmallWorld,latticeDim = 1,numberOfEdgesDelRandomly =numberOFDelforSmallWorld,p1=Probability,p2=Probability,PstepSize =1 ))
    print(paste("Serializing and saving Network objects Object:SmallWorldEdges",Probability,sep = ""))
    saveRDS(SmallWorldEdges, file = paste(savingDir,"SmallWorld",Probability,sep = ""), ascii = FALSE, version = NULL,
            compress = TRUE, refhook = NULL)

  }

  #############
  if(plotGraph=="Y")
  {
   # require('ggplot2')
    ########Plotting
    DfGlobalEdgeVertices<-NULL

    DfEdgeGlobalClusteringCoefficent<-NULL
    DfGlobalEdgeCentralityBetweenessMean<-NULL
    DfGlobalEdgeCentralityClosenessMean<-NULL
    DfGlobalEdgeCentralityDegreeMean<-NULL
    DfGlobalEdgeAvgGeodesicPath<-NULL




    print("loading serialized object-randomGraphEdgeComparisn from hdd")
    randomGraphEdgeComparisn = readRDS(paste(savingDir,"randomGraphEdgeComparisn",sep =""), refhook = NULL)

    dfRandomTemp = Plot2dListOfRandomGraphPropertiesMean(randomGraphEdgeComparisn,x="NumberOfEdges",y="GlobalClusteringCoefficent",xlabel = "Number of Edges",ylabel = "Global Clustering Coefficent")
    dfRandomTemp = cbind(dfRandomTemp,networkType="randomGraph")
    DfEdgeGlobalClusteringCoefficent<-rbind(DfEdgeGlobalClusteringCoefficent,dfRandomTemp)
    dfRandomTemp=NULL

    dfRandomTemp = Plot2dListOfRandomGraphPropertiesMean(randomGraphEdgeComparisn,x="NumberOfEdges",y="Vertices",xlabel = "Number of Edges",ylabel = "Vertices")
    dfRandomTemp = cbind(dfRandomTemp,networkType="randomGraph")
    DfGlobalEdgeVertices<-rbind(DfGlobalEdgeVertices,dfRandomTemp)
    dfRandomTemp=NULL

    dfRandomTemp = Plot2dListOfRandomGraphPropertiesMean(randomGraphEdgeComparisn,x="NumberOfEdges",y="CentralityBetweenessMean",xlabel = "Number of Edges",ylabel = "Mean Betweeness Centrality")
    dfRandomTemp = cbind(dfRandomTemp,networkType="randomGraph")
    DfGlobalEdgeCentralityBetweenessMean<-rbind(DfGlobalEdgeCentralityBetweenessMean,dfRandomTemp)
    dfRandomTemp=NULL

    dfRandomTemp = Plot2dListOfRandomGraphPropertiesMean(randomGraphEdgeComparisn,x="NumberOfEdges",y="CentralityClosenessMean",xlabel = "Number of Edges",ylabel = "Mean Closeness Centrality")
    dfRandomTemp = cbind(dfRandomTemp,networkType="randomGraph")
    DfGlobalEdgeCentralityClosenessMean<-rbind(DfGlobalEdgeCentralityClosenessMean,dfRandomTemp)
    dfRandomTemp=NULL

    dfRandomTemp = Plot2dListOfRandomGraphPropertiesMean(randomGraphEdgeComparisn,x="NumberOfEdges",y="CentralityDegreeMean",xlabel = "Number of Edges",ylabel = "Mean Degree Centrality")
    dfRandomTemp = cbind(dfRandomTemp,networkType="randomGraph")
    DfGlobalEdgeCentralityDegreeMean<-rbind(DfGlobalEdgeCentralityDegreeMean,dfRandomTemp)
    dfRandomTemp=NULL


    dfRandomTemp = Plot2dListOfRandomGraphPropertiesMean(randomGraphEdgeComparisn,x="NumberOfEdges",y="AvgGeodesicPath",xlabel = "Number of Edges",ylabel = "Avg Geodesic Path")
    dfRandomTemp = cbind(dfRandomTemp,networkType="randomGraph")
    DfGlobalEdgeAvgGeodesicPath<-rbind(DfGlobalEdgeAvgGeodesicPath,dfRandomTemp)
    dfRandomTemp=NULL







    for (powersf in ScaleFreePowerRange) {
      for (probabilitysw in SmallWorldProbability) {


        NameScaleFree <- paste("Scalefree",powersf,sep ="")
        NameSmallWorld <- paste("SmallWorld",probabilitysw,sep ="")
        NameFileName <- paste("sc",powersf,"sw",probabilitysw,sep ="")

        print(paste("loading serialized object from hdd:",NameScaleFree))
        Scalefree = readRDS(paste(savingDir,NameScaleFree,sep =""), refhook = NULL)

        print(paste("loading serialized object from hdd:",NameSmallWorld))
        SmallWorldEdges = readRDS(paste(savingDir,NameSmallWorld,sep =""), refhook = NULL)


        dfScaleFreeMultiPlot = Plot2dListOfRandomGraphPropertiesMean(Scalefree,x="NumberofEdges",y="GlobalClusteringCoefficent",xlabel = "Number of Edges",ylabel = "Global Clustering Coefficent")
        dfScaleFreeMultiPlot = cbind(dfScaleFreeMultiPlot,networkType=NameScaleFree)
        DfEdgeGlobalClusteringCoefficent<-rbind(DfEdgeGlobalClusteringCoefficent,dfScaleFreeMultiPlot)



        dfSmallWorldMultiPlot = Plot2dListOfRandomGraphPropertiesMean(SmallWorldEdges,x="NumberOfEdges",y="GlobalClusteringCoefficent",xlabel = "Number of Edges",ylabel = "Global Clustering Coefficent")
        dfSmallWorldMultiPlot = cbind(dfSmallWorldMultiPlot,networkType=NameSmallWorld)
        DfEdgeGlobalClusteringCoefficent<-rbind(DfEdgeGlobalClusteringCoefficent,dfSmallWorldMultiPlot)


        dfScaleFreeMultiPlot = Plot2dListOfRandomGraphPropertiesMean(Scalefree,x="NumberofEdges",y="Vertices",xlabel = "Number of Edges",ylabel = "Vertices")
        dfScaleFreeMultiPlot = cbind(dfScaleFreeMultiPlot,networkType=NameScaleFree)
        DfGlobalEdgeVertices<-rbind(DfGlobalEdgeVertices,dfScaleFreeMultiPlot)



        dfSmallWorldMultiPlot = Plot2dListOfRandomGraphPropertiesMean(SmallWorldEdges,x="NumberOfEdges",y="Vertices",xlabel = "Number of Edges",ylabel = "Vertices")
        dfSmallWorldMultiPlot = cbind(dfSmallWorldMultiPlot,networkType=NameSmallWorld)
        DfGlobalEdgeVertices<-rbind(DfGlobalEdgeVertices,dfSmallWorldMultiPlot)



        #print(paste("generating and saving plots for:",NameFileName))

        dfRandomTemp = Plot2dListOfRandomGraphPropertiesMean(Scalefree,x="NumberofEdges",y="CentralityBetweenessMean",xlabel = "Number of Edges",ylabel = "Mean Betweeness Centrality")
        dfRandomTemp = cbind(dfRandomTemp,networkType=NameScaleFree)
        DfGlobalEdgeCentralityBetweenessMean<-rbind(DfGlobalEdgeCentralityBetweenessMean,dfRandomTemp)
        dfRandomTemp=NULL

        dfRandomTemp = Plot2dListOfRandomGraphPropertiesMean(Scalefree,x="NumberofEdges",y="CentralityClosenessMean",xlabel = "Number of Edges",ylabel = "Mean Closeness Centrality")
        dfRandomTemp = cbind(dfRandomTemp,networkType=NameScaleFree)
        DfGlobalEdgeCentralityClosenessMean<-rbind(DfGlobalEdgeCentralityClosenessMean,dfRandomTemp)
        dfRandomTemp=NULL

        dfRandomTemp = Plot2dListOfRandomGraphPropertiesMean(Scalefree,x="NumberofEdges",y="CentralityDegreeMean",xlabel = "Number of Edges",ylabel = "Mean Degree Centrality")
        dfRandomTemp = cbind(dfRandomTemp,networkType=NameScaleFree)
        DfGlobalEdgeCentralityDegreeMean<-rbind(DfGlobalEdgeCentralityDegreeMean,dfRandomTemp)
        dfRandomTemp=NULL


        dfRandomTemp = Plot2dListOfRandomGraphPropertiesMean(Scalefree,x="NumberofEdges",y="AvgGeodesicPath",xlabel = "Number of Edges",ylabel = "Avg Geodesic Path")
        dfRandomTemp = cbind(dfRandomTemp,networkType=NameScaleFree)
        DfGlobalEdgeAvgGeodesicPath<-rbind(DfGlobalEdgeAvgGeodesicPath,dfRandomTemp)
        dfRandomTemp=NULL


        #small world

        dfRandomTemp = Plot2dListOfRandomGraphPropertiesMean(SmallWorldEdges,x="NumberOfEdges",y="CentralityBetweenessMean",xlabel = "Number of Edges",ylabel = "Mean Betweeness Centrality")
        dfRandomTemp = cbind(dfRandomTemp,networkType=NameSmallWorld)
        DfGlobalEdgeCentralityBetweenessMean<-rbind(DfGlobalEdgeCentralityBetweenessMean,dfRandomTemp)
        dfRandomTemp=NULL

        dfRandomTemp = Plot2dListOfRandomGraphPropertiesMean(SmallWorldEdges,x="NumberOfEdges",y="CentralityClosenessMean",xlabel = "Number of Edges",ylabel = "Mean Closeness Centrality")
        dfRandomTemp = cbind(dfRandomTemp,networkType=NameSmallWorld)
        DfGlobalEdgeCentralityClosenessMean<-rbind(DfGlobalEdgeCentralityClosenessMean,dfRandomTemp)
        dfRandomTemp=NULL

        dfRandomTemp = Plot2dListOfRandomGraphPropertiesMean(SmallWorldEdges,x="NumberOfEdges",y="CentralityDegreeMean",xlabel = "Number of Edges",ylabel = "Mean Degree Centrality")
        dfRandomTemp = cbind(dfRandomTemp,networkType=NameSmallWorld)
        DfGlobalEdgeCentralityDegreeMean<-rbind(DfGlobalEdgeCentralityDegreeMean,dfRandomTemp)
        dfRandomTemp=NULL


        dfRandomTemp = Plot2dListOfRandomGraphPropertiesMean(SmallWorldEdges,x="NumberOfEdges",y="AvgGeodesicPath",xlabel = "Number of Edges",ylabel = "Avg Geodesic Path")
        dfRandomTemp = cbind(dfRandomTemp,networkType=NameSmallWorld)
        DfGlobalEdgeAvgGeodesicPath<-rbind(DfGlobalEdgeAvgGeodesicPath,dfRandomTemp)
        dfRandomTemp=NULL






      }


    }


    p1 <- ggplot(DfEdgeGlobalClusteringCoefficent, aes(`Number of Edges`,`Global Clustering Coefficent`, group = networkType,
                                                       colour = networkType)) + geom_line(size = 1)
    #p1 + geom_text(data = DfEdgeGlobalClusteringCoefficent, aes(label = networkType), hjust = 0.7, vjust = 1)

    ggsave(paste(savingDir,"EdgeGlobalClusteringCoef.pdf",sep=""), width = 20, height = 20, units = "cm")


    p1 <- ggplot(DfGlobalEdgeVertices, aes(`Number of Edges`,Vertices, group = networkType,
                                           colour = networkType)) + geom_line(size = 1)
    #p1 + geom_text(data = DfGlobalEdgeVertices, aes(label = networkType), hjust = 0.7, vjust = 1)

    ggsave(paste(savingDir,"EdgeVertices.pdf"), width = 20, height = 20, units = "cm")


    p1 <- ggplot(DfGlobalEdgeCentralityDegreeMean, aes(`Number of Edges`,`Mean Degree Centrality`, group = networkType,
                                                       colour = networkType)) + geom_line(size = 1)
    #p1 + geom_text(data = DfGlobalEdgeCentralityDegreesMean, aes(label = networkType), hjust = 0.7, vjust = 1)

    ggsave(paste(savingDir,"EdgeMeanDegreeCentrality.pdf",sep=""), width = 20, height = 20, units = "cm")


    p1 <- ggplot(DfGlobalEdgeCentralityBetweenessMean, aes(`Number of Edges`,`Mean Betweeness Centrality`, group = networkType,
                                                           colour = networkType)) + geom_line(size = 1)
    #p1 + geom_text(data = DfGlobalEdgeCentralityBetweenessMean, aes(label = networkType), hjust = 0.7, vjust = 1)

    ggsave(paste(savingDir,"EdgeMeanBetweenessCentrality.pdf",sep=""), width = 20, height = 20, units = "cm")

    p1 <- ggplot(DfGlobalEdgeCentralityClosenessMean, aes(`Number of Edges`,`Mean Closeness Centrality`, group = networkType,
                                                          colour = networkType)) + geom_line(size = 1)
    #p1 + geom_text(data = DfGlobalEdgeClosenessMean, aes(label = networkType), hjust = 0.7, vjust = 1)

    ggsave(paste(savingDir,"EdgeMeanClosnessCentrality.pdf",sep=""), width = 20, height = 20, units = "cm")


    p1 <- ggplot(DfGlobalEdgeAvgGeodesicPath, aes(`Number of Edges`,`Avg Geodesic Path`, group = networkType,
                                                  colour = networkType)) + geom_line(size = 1)
    #p1 + geom_text(data = DfGlobalEdgeAvgGeodesicPath, aes(label = networkType), hjust = 0.7, vjust = 1)

    ggsave(paste(savingDir,"EdgeAvgGeodesicPath.pdf",sep=""), width = 20, height = 20, units = "cm")


  }


  print("Completed!( Number of Edges = Number of Vertices*2-3, this is to match the number of vertices in all three graps) ")
  # Stop the clock
  print("user time:execution of the code, system time:CPU ,elapsed time: total")
  proc.time() - ptm


}

#'
#'
#' This function allowes you to generate three different types of networks with same number of vertices and edges
#' @param ScaleFreePowerRange A vector containing powers for scale free networks, e.g. c(1.5,1.75)
#' @param SmallWorldProbability A vector containing probabilities for small world networks, e.g. c(0.5,0.6)
#' @param plot type of plot(string input "point" or "line")
#' @param savingDir A string with the directory where to save graph plots. This function requires the network proporties to be generated and saved in that folder by using gen_graps function
#' @import igraph
#' @import ggplot2
#' @export
gen_plots<- function(ScaleFreePowerRange,SmallWorldProbability,plot=C("point","line"),savingDir)
{
  #  require('igraph')
  #  require('ggplot2')
  # Start the clock!
  ptm <- proc.time()

  print("generating Scalefree networks(preferential attachment, 2 edges each step) and calculating different properties")

  # EdgesVector = (VerticesVector * 2)-3

  #EdgesVector = ((VerticesVector*edgesEachScaleFree)-((edgesEachScaleFree*(edgesEachScaleFree+1))/2))

  #  NeiSmallWorld = edgesEachScaleFree

  #  numberOFDelforSmallWorld = (VerticesVector*(edgesEachScaleFree-NeiSmallWorld))+((edgesEachScaleFree*(edgesEachScaleFree+1))/2)
  # print(numberOFDelforSmallWorld)

  ##numberOFDelforSmallWorld = ((edgesEachScaleFree*(edgesEachScaleFree+1))/2)


  #############

  ########Plotting
  DfGlobalEdgeVertices<-NULL

  DfEdgeGlobalClusteringCoefficent<-NULL
  DfGlobalEdgeCentralityBetweenessMean<-NULL
  DfGlobalEdgeCentralityClosenessMean<-NULL
  DfGlobalEdgeCentralityDegreeMean<-NULL
  DfGlobalEdgeAvgGeodesicPath<-NULL




  print("loading serialized object-randomGraphEdgeComparisn from hdd")
  randomGraphEdgeComparisn = readRDS(paste(savingDir,"randomGraphEdgeComparisn",sep =""), refhook = NULL)

  dfRandomTemp = Plot2dListOfRandomGraphPropertiesMean(randomGraphEdgeComparisn,x="NumberOfEdges",y="GlobalClusteringCoefficent",xlabel = "Number of Edges",ylabel = "Global Clustering Coefficent")
  dfRandomTemp = cbind(dfRandomTemp,networkType="randomGraph")
  DfEdgeGlobalClusteringCoefficent<-rbind(DfEdgeGlobalClusteringCoefficent,dfRandomTemp)
  dfRandomTemp=NULL

  dfRandomTemp = Plot2dListOfRandomGraphPropertiesMean(randomGraphEdgeComparisn,x="NumberOfEdges",y="Vertices",xlabel = "Number of Edges",ylabel = "Vertices")
  dfRandomTemp = cbind(dfRandomTemp,networkType="randomGraph")
  DfGlobalEdgeVertices<-rbind(DfGlobalEdgeVertices,dfRandomTemp)
  dfRandomTemp=NULL

  dfRandomTemp = Plot2dListOfRandomGraphPropertiesMean(randomGraphEdgeComparisn,x="NumberOfEdges",y="CentralityBetweenessMean",xlabel = "Number of Edges",ylabel = "Mean Betweeness Centrality")
  dfRandomTemp = cbind(dfRandomTemp,networkType="randomGraph")
  DfGlobalEdgeCentralityBetweenessMean<-rbind(DfGlobalEdgeCentralityBetweenessMean,dfRandomTemp)
  dfRandomTemp=NULL

  dfRandomTemp = Plot2dListOfRandomGraphPropertiesMean(randomGraphEdgeComparisn,x="NumberOfEdges",y="CentralityClosenessMean",xlabel = "Number of Edges",ylabel = "Mean Closeness Centrality")
  dfRandomTemp = cbind(dfRandomTemp,networkType="randomGraph")
  DfGlobalEdgeCentralityClosenessMean<-rbind(DfGlobalEdgeCentralityClosenessMean,dfRandomTemp)
  dfRandomTemp=NULL

  dfRandomTemp = Plot2dListOfRandomGraphPropertiesMean(randomGraphEdgeComparisn,x="NumberOfEdges",y="CentralityDegreeMean",xlabel = "Number of Edges",ylabel = "Mean Degree Centrality")
  dfRandomTemp = cbind(dfRandomTemp,networkType="randomGraph")
  DfGlobalEdgeCentralityDegreeMean<-rbind(DfGlobalEdgeCentralityDegreeMean,dfRandomTemp)
  dfRandomTemp=NULL


  dfRandomTemp = Plot2dListOfRandomGraphPropertiesMean(randomGraphEdgeComparisn,x="NumberOfEdges",y="AvgGeodesicPath",xlabel = "Number of Edges",ylabel = "Avg Geodesic Path")
  dfRandomTemp = cbind(dfRandomTemp,networkType="randomGraph")
  DfGlobalEdgeAvgGeodesicPath<-rbind(DfGlobalEdgeAvgGeodesicPath,dfRandomTemp)
  dfRandomTemp=NULL







  for (powersf in ScaleFreePowerRange) {
    for (probabilitysw in SmallWorldProbability) {


      NameScaleFree <- paste("Scalefree",powersf,sep ="")
      NameSmallWorld <- paste("SmallWorld",probabilitysw,sep ="")
      NameFileName <- paste("sc",powersf,"sw",probabilitysw,sep ="")

      print(paste("loading serialized object from hdd:",NameScaleFree))
      Scalefree = readRDS(paste(savingDir,NameScaleFree,sep =""), refhook = NULL)

      print(paste("loading serialized object from hdd:",NameSmallWorld))
      SmallWorldEdges = readRDS(paste(savingDir,NameSmallWorld,sep =""), refhook = NULL)


      dfScaleFreeMultiPlot = Plot2dListOfRandomGraphPropertiesMean(Scalefree,x="NumberofEdges",y="GlobalClusteringCoefficent",xlabel = "Number of Edges",ylabel = "Global Clustering Coefficent")
      dfScaleFreeMultiPlot = cbind(dfScaleFreeMultiPlot,networkType=NameScaleFree)
      DfEdgeGlobalClusteringCoefficent<-rbind(DfEdgeGlobalClusteringCoefficent,dfScaleFreeMultiPlot)



      dfSmallWorldMultiPlot = Plot2dListOfRandomGraphPropertiesMean(SmallWorldEdges,x="NumberOfEdges",y="GlobalClusteringCoefficent",xlabel = "Number of Edges",ylabel = "Global Clustering Coefficent")
      dfSmallWorldMultiPlot = cbind(dfSmallWorldMultiPlot,networkType=NameSmallWorld)
      DfEdgeGlobalClusteringCoefficent<-rbind(DfEdgeGlobalClusteringCoefficent,dfSmallWorldMultiPlot)


      dfScaleFreeMultiPlot = Plot2dListOfRandomGraphPropertiesMean(Scalefree,x="NumberofEdges",y="Vertices",xlabel = "Number of Edges",ylabel = "Vertices")
      dfScaleFreeMultiPlot = cbind(dfScaleFreeMultiPlot,networkType=NameScaleFree)
      DfGlobalEdgeVertices<-rbind(DfGlobalEdgeVertices,dfScaleFreeMultiPlot)



      dfSmallWorldMultiPlot = Plot2dListOfRandomGraphPropertiesMean(SmallWorldEdges,x="NumberOfEdges",y="Vertices",xlabel = "Number of Edges",ylabel = "Vertices")
      dfSmallWorldMultiPlot = cbind(dfSmallWorldMultiPlot,networkType=NameSmallWorld)
      DfGlobalEdgeVertices<-rbind(DfGlobalEdgeVertices,dfSmallWorldMultiPlot)



      #print(paste("generating and saving plots for:",NameFileName))

      dfRandomTemp = Plot2dListOfRandomGraphPropertiesMean(Scalefree,x="NumberofEdges",y="CentralityBetweenessMean",xlabel = "Number of Edges",ylabel = "Mean Betweeness Centrality")
      dfRandomTemp = cbind(dfRandomTemp,networkType=NameScaleFree)
      DfGlobalEdgeCentralityBetweenessMean<-rbind(DfGlobalEdgeCentralityBetweenessMean,dfRandomTemp)
      dfRandomTemp=NULL

      dfRandomTemp = Plot2dListOfRandomGraphPropertiesMean(Scalefree,x="NumberofEdges",y="CentralityClosenessMean",xlabel = "Number of Edges",ylabel = "Mean Closeness Centrality")
      dfRandomTemp = cbind(dfRandomTemp,networkType=NameScaleFree)
      DfGlobalEdgeCentralityClosenessMean<-rbind(DfGlobalEdgeCentralityClosenessMean,dfRandomTemp)
      dfRandomTemp=NULL

      dfRandomTemp = Plot2dListOfRandomGraphPropertiesMean(Scalefree,x="NumberofEdges",y="CentralityDegreeMean",xlabel = "Number of Edges",ylabel = "Mean Degree Centrality")
      dfRandomTemp = cbind(dfRandomTemp,networkType=NameScaleFree)
      DfGlobalEdgeCentralityDegreeMean<-rbind(DfGlobalEdgeCentralityDegreeMean,dfRandomTemp)
      dfRandomTemp=NULL


      dfRandomTemp = Plot2dListOfRandomGraphPropertiesMean(Scalefree,x="NumberofEdges",y="AvgGeodesicPath",xlabel = "Number of Edges",ylabel = "Avg Geodesic Path")
      dfRandomTemp = cbind(dfRandomTemp,networkType=NameScaleFree)
      DfGlobalEdgeAvgGeodesicPath<-rbind(DfGlobalEdgeAvgGeodesicPath,dfRandomTemp)
      dfRandomTemp=NULL


      #small world

      dfRandomTemp = Plot2dListOfRandomGraphPropertiesMean(SmallWorldEdges,x="NumberOfEdges",y="CentralityBetweenessMean",xlabel = "Number of Edges",ylabel = "Mean Betweeness Centrality")
      dfRandomTemp = cbind(dfRandomTemp,networkType=NameSmallWorld)
      DfGlobalEdgeCentralityBetweenessMean<-rbind(DfGlobalEdgeCentralityBetweenessMean,dfRandomTemp)
      dfRandomTemp=NULL

      dfRandomTemp = Plot2dListOfRandomGraphPropertiesMean(SmallWorldEdges,x="NumberOfEdges",y="CentralityClosenessMean",xlabel = "Number of Edges",ylabel = "Mean Closeness Centrality")
      dfRandomTemp = cbind(dfRandomTemp,networkType=NameSmallWorld)
      DfGlobalEdgeCentralityClosenessMean<-rbind(DfGlobalEdgeCentralityClosenessMean,dfRandomTemp)
      dfRandomTemp=NULL

      dfRandomTemp = Plot2dListOfRandomGraphPropertiesMean(SmallWorldEdges,x="NumberOfEdges",y="CentralityDegreeMean",xlabel = "Number of Edges",ylabel = "Mean Degree Centrality")
      dfRandomTemp = cbind(dfRandomTemp,networkType=NameSmallWorld)
      DfGlobalEdgeCentralityDegreeMean<-rbind(DfGlobalEdgeCentralityDegreeMean,dfRandomTemp)
      dfRandomTemp=NULL


      dfRandomTemp = Plot2dListOfRandomGraphPropertiesMean(SmallWorldEdges,x="NumberOfEdges",y="AvgGeodesicPath",xlabel = "Number of Edges",ylabel = "Avg Geodesic Path")
      dfRandomTemp = cbind(dfRandomTemp,networkType=NameSmallWorld)
      DfGlobalEdgeAvgGeodesicPath<-rbind(DfGlobalEdgeAvgGeodesicPath,dfRandomTemp)
      dfRandomTemp=NULL






    }

  }


  if(plot=="line")
  {
    p1 <- ggplot(DfEdgeGlobalClusteringCoefficent, aes(`Number of Edges`,`Global Clustering Coefficent`, group = networkType,
                                                       colour = networkType)) + geom_line(size = 1)
    #p1 + geom_text(data = DfEdgeGlobalClusteringCoefficent, aes(label = networkType), hjust = 0.7, vjust = 1)

    ggsave(paste(savingDir,"EdgeGlobalClusteringCoef.pdf",sep=""), width = 20, height = 20, units = "cm")


    p1 <- ggplot(DfGlobalEdgeVertices, aes(`Number of Edges`,Vertices, group = networkType,
                                           colour = networkType)) + geom_line(size = 1)
    #p1 + geom_text(data = DfGlobalEdgeVertices, aes(label = networkType), hjust = 0.7, vjust = 1)

    ggsave(paste(savingDir,"EdgeVertices.pdf"), width = 20, height = 20, units = "cm")


    p1 <- ggplot(DfGlobalEdgeCentralityDegreeMean, aes(`Number of Edges`,`Mean Degree Centrality`, group = networkType,
                                                       colour = networkType)) + geom_line(size = 1)
    #p1 + geom_text(data = DfGlobalEdgeCentralityDegreesMean, aes(label = networkType), hjust = 0.7, vjust = 1)

    ggsave(paste(savingDir,"EdgeMeanDegreeCentrality.pdf",sep=""), width = 20, height = 20, units = "cm")


    p1 <- ggplot(DfGlobalEdgeCentralityBetweenessMean, aes(`Number of Edges`,`Mean Betweeness Centrality`, group = networkType,
                                                           colour = networkType)) + geom_line(size = 1)
    #p1 + geom_text(data = DfGlobalEdgeCentralityBetweenessMean, aes(label = networkType), hjust = 0.7, vjust = 1)

    ggsave(paste(savingDir,"EdgeMeanBetweenessCentrality.pdf",sep=""), width = 20, height = 20, units = "cm")

    p1 <- ggplot(DfGlobalEdgeCentralityClosenessMean, aes(`Number of Edges`,`Mean Closeness Centrality`, group = networkType,
                                                          colour = networkType)) + geom_line(size = 1)
    #p1 + geom_text(data = DfGlobalEdgeClosenessMean, aes(label = networkType), hjust = 0.7, vjust = 1)

    ggsave(paste(savingDir,"EdgeMeanClosnessCentrality.pdf",sep=""), width = 20, height = 20, units = "cm")


    p1 <- ggplot(DfGlobalEdgeAvgGeodesicPath, aes(`Number of Edges`,`Avg Geodesic Path`, group = networkType,
                                                  colour = networkType)) + geom_line(size = 1)
    #p1 + geom_text(data = DfGlobalEdgeAvgGeodesicPath, aes(label = networkType), hjust = 0.7, vjust = 1)

    ggsave(paste(savingDir,"EdgeAvgGeodesicPath.pdf",sep=""), width = 20, height = 20, units = "cm")

  }
  if(plot=="point"){

    p1 <- ggplot(DfEdgeGlobalClusteringCoefficent, aes(`Number of Edges`,`Global Clustering Coefficent`, group = networkType,
                                                       colour = networkType)) + geom_point(size = 1)
    #p1 + geom_text(data = DfEdgeGlobalClusteringCoefficent, aes(label = networkType), hjust = 0.7, vjust = 1)

    ggsave(paste(savingDir,plot,"EdgeGlobalClusteringCoef.pdf",sep=""), width = 20, height = 20, units = "cm")


    p1 <- ggplot(DfGlobalEdgeVertices, aes(`Number of Edges`,Vertices, group = networkType,
                                           colour = networkType)) + geom_point(size = 1)
    #p1 + geom_text(data = DfGlobalEdgeVertices, aes(label = networkType), hjust = 0.7, vjust = 1)

    ggsave(paste(savingDir,plot,"EdgeVertices.pdf"), width = 20, height = 20, units = "cm")


    p1 <- ggplot(DfGlobalEdgeCentralityDegreeMean, aes(`Number of Edges`,`Mean Degree Centrality`, group = networkType,
                                                       colour = networkType)) + geom_point(size = 1)
    #p1 + geom_text(data = DfGlobalEdgeCentralityDegreesMean, aes(label = networkType), hjust = 0.7, vjust = 1)

    ggsave(paste(savingDir,plot,"EdgeMeanDegreeCentrality.pdf",sep=""), width = 20, height = 20, units = "cm")


    p1 <- ggplot(DfGlobalEdgeCentralityBetweenessMean, aes(`Number of Edges`,`Mean Betweeness Centrality`, group = networkType,
                                                           colour = networkType)) + geom_point(size = 1)
    #p1 + geom_text(data = DfGlobalEdgeCentralityBetweenessMean, aes(label = networkType), hjust = 0.7, vjust = 1)

    ggsave(paste(savingDir,plot,"EdgeMeanBetweenessCentrality.pdf",sep=""), width = 20, height = 20, units = "cm")

    p1 <- ggplot(DfGlobalEdgeCentralityClosenessMean, aes(`Number of Edges`,`Mean Closeness Centrality`, group = networkType,
                                                          colour = networkType)) + geom_point(size = 1)
    #p1 + geom_text(data = DfGlobalEdgeClosenessMean, aes(label = networkType), hjust = 0.7, vjust = 1)

    ggsave(paste(savingDir,plot,"EdgeMeanClosnessCentrality.pdf",sep=""), width = 20, height = 20, units = "cm")


    p1 <- ggplot(DfGlobalEdgeAvgGeodesicPath, aes(`Number of Edges`,`Avg Geodesic Path`, group = networkType,
                                                  colour = networkType)) + geom_point(size = 1)
    #p1 + geom_text(data = DfGlobalEdgeAvgGeodesicPath, aes(label = networkType), hjust = 0.7, vjust = 1)

    ggsave(paste(savingDir,plot,"EdgeAvgGeodesicPath.pdf",sep=""), width = 20, height = 20, units = "cm")

  }



  print("Completed!( Number of Edges = Number of Vertices*2-3, this is to match the number of vertices in all three graps) ")
  # Stop the clock
  print("user time:execution of the code, system time:CPU ,elapsed time: total")
  proc.time() - ptm


}



createDataFrameforPlotting<-function(dfObj,name)
{

  dfObj = cbind(dfObj,graphType=name)




}





DelRandomEdge <- function(graphObj)
{

  gsize(graphObj)
  edgeList = get.edgelist(graphObj)
  randomNumberlength=length(edgeList)/2
  rs <- sample_seq(1, randomNumberlength, 1)
  SelectedEdge = edgeList[rs,]
  EdgeforDel = paste((SelectedEdge)[1],"|",(SelectedEdge)[2],sep = "")
  graphObj = graphObj%>%delete_edges(EdgeforDel)
  gsize(graphObj)


  return(graphObj)

}


ERGraphRandomPropertiesGNM2<- function(n1=NULL,p1=NULL,n2=NULL,p2=NULL,VectorVertices=NULL,VectorEdges=NULL,PstepSize,sampleSize)
{
  #	Measure other properties of the graphs:
  #- Average geodesic/shortest path, global clustering coefficient, degree centrality(degree, closeness, betweeness)

  Probability<-1
  Vertices<-1
  NumberOfEdges<-1

  CentralityDegreeMean<- 1
  CentralityClosenessMean<-1
  CentralityBetweenessMean<- 1
  CentralityDegreeMin<- 1
  CentralityClosenessMin<-1
  CentralityBetweenessmin<- 1
  CentralityDegreeMax<- 1
  CentralityClosenessmax<-1
  CentralityBetweenessmax<- 1

  CentralityDegreeList<-1
  names(CentralityDegreeList) = "Degree"
  CentralityClosenessList<-1
  names(CentralityClosenessList) = "Closeness"
  CentralityBetweenessList<-1
  names(CentralityBetweenessList) = "Estimate_betweenness"

  GlobalClusteringCoefficent <-1
  AvgGeodesicPath<- 1
  NumberOfEdges<- 1

  counter1 = 1

  for (n in VectorVertices)
  {
    e=VectorEdges[[counter1]]

    erRandGraph = erdos.renyi.game(n=n,p=e, type = "gnm", directed = FALSE,loops = FALSE)
    ##print(paste(" Edges ",e))
    # #print(paste(" vertices ",n))
    #calculation of properties from graph

    Degree = degree(erRandGraph, normalized = FALSE)
    Closeness  = closeness.estimate(erRandGraph, cutoff = -1)
    Estimate_betweenness  = estimate_betweenness(erRandGraph, cutoff = -1)
    Mean_distance = mean_distance(erRandGraph, directed = FALSE)
    Transitivity = transitivity(erRandGraph, type = c("global"), vids = NULL, weights = NULL)

    # Post calculation of properties

    CentralityDegreeList[counter1]=list(Degree)

    CentralityClosenessList[counter1] = list(Closeness)

    CentralityBetweenessList[counter1] = list(Estimate_betweenness)

    CentralityDegreeMean[[counter1]] =  mean(Degree)
    CentralityClosenessMean[[counter1]] =  mean(Closeness)
    CentralityBetweenessMean[[counter1]] =  mean(Estimate_betweenness) #cutoff:The maximum path length to consider when calculating the betweenness. If zero or negative then there is no such limit.
    CentralityDegreeMin[[counter1]] =  min(Degree)
    CentralityClosenessMin[[counter1]] =  min(Closeness)
    CentralityBetweenessmin[[counter1]] =  min(Estimate_betweenness)
    CentralityDegreeMax[[counter1]] =  max(Degree)
    CentralityClosenessmax[[counter1]] =  max(Closeness)
    CentralityBetweenessmax[[counter1]] =  max(Estimate_betweenness)

    AvgGeodesicPath[[counter1]] = Mean_distance
    GlobalClusteringCoefficent[[counter1]] = Transitivity

    #Graph Info

    Vertices[[counter1]]<-n
    NumberOfEdges[[counter1]]<-e

    #Next graph counter increament

    counter1 = counter1+1;


  }




  #return(CentralityDegreeList)

  return(list(data.frame(Vertices,NumberOfEdges,GlobalClusteringCoefficent,CentralityDegreeMean,CentralityClosenessMean,CentralityBetweenessMean,AvgGeodesicPath),CentralityDegreeList,CentralityClosenessList));
  #return(data.frame(Vertices,Probability,GlobalClusteringCoefficent,CentralityDegreeMean,CentralityClosenessMean,CentralityBetweenessMean,AvgGeodesicPath, check.rows = FALSE));



}



BAGraphScalefreePropertiesEdges<- function(n1=NULL,p1,n2=NULL,p2,PstepSize,VectorVertices=NULL, sampleSize,EdgesEachStep=1)
{
  #	Measure other properties of the graphs:
  #- Average geodesic/shortest path, global clustering coefficient, degree centrality(degree, closeness, betweeness)

  Probability<-1
  Vertices<-1

  CentralityDegreeMean<- 1
  CentralityClosenessMean<-1
  CentralityBetweenessMean<- 1
  CentralityDegreeMin<- 1
  CentralityClosenessMin<-1
  CentralityBetweenessmin<- 1
  CentralityDegreeMax<- 1
  CentralityClosenessmax<-1
  CentralityBetweenessmax<- 1

  CentralityDegreeList<-1
  names(CentralityDegreeList) = "Degree"
  CentralityClosenessList<-1
  names(CentralityClosenessList) = "Closeness"
  CentralityBetweenessList<-1
  names(CentralityBetweenessList) = "Estimate_betweenness"

  GlobalClusteringCoefficent <-1
  AvgGeodesicPath<- 1
  NumberofEdges<- 1

  counter1 = 1

  if(is.null(VectorVertices)){
    for (n in n1:n2)
    {
      flag = 0;
      p=p1;

      repeat {

        if(p>p2)
        {
          break
        }
        Scalefree = sample_pa(n, power = p, m = EdgesEachStep, out.dist = NULL, out.seq = NULL,
                              out.pref = FALSE, zero.appeal = 1, directed = FALSE,
                              algorithm = "psumtree", start.graph = NULL)
        ##print(paste(" probablity ",p))
        # #print(paste(" vertices ",n))

        #calculation of properties from graph

        Degree = degree(Scalefree, normalized = FALSE)
        Closeness  = closeness.estimate(Scalefree, cutoff = -1)
        Estimate_betweenness  = estimate_betweenness(Scalefree, cutoff = -1)
        Mean_distance = mean_distance(Scalefree, directed = FALSE)
        Transitivity = transitivity(Scalefree, type = c("global"), vids = NULL, weights = NULL)
        Edges_get= gsize(Scalefree)
        # Post calculation of properties
        #print(paste(" Edges:",Edges_get))
        CentralityDegreeList[counter1]=list(Degree)

        CentralityClosenessList[counter1] = list(Closeness)

        CentralityBetweenessList[counter1] = list(Estimate_betweenness)

        CentralityDegreeMean[[counter1]] =  mean(Degree)
        CentralityClosenessMean[[counter1]] =  mean(Closeness)
        CentralityBetweenessMean[[counter1]] =  mean(Estimate_betweenness) #cutoff:The maximum path length to consider when calculating the betweenness. If zero or negative then there is no such limit.
        CentralityDegreeMin[[counter1]] =  min(Degree)
        CentralityClosenessMin[[counter1]] =  min(Closeness)
        CentralityBetweenessmin[[counter1]] =  min(Estimate_betweenness)
        CentralityDegreeMax[[counter1]] =  max(Degree)
        CentralityClosenessmax[[counter1]] =  max(Closeness)
        CentralityBetweenessmax[[counter1]] =  max(Estimate_betweenness)

        AvgGeodesicPath[[counter1]] = Mean_distance
        GlobalClusteringCoefficent[[counter1]] = Transitivity

        #Graph Info

        Vertices[[counter1]]<-n
        Probability[[counter1]]<-p
        NumberofEdges[[counter1]]<- Edges_get

        #Next graph counter increament

        counter1 = counter1+1;
        p= p+PstepSize;

      }

    }
  }
  else{
    for (n in VectorVertices)
    {
      {
        flag = 0;
        p=p1;

        repeat {

          if(p>p2)
          {
            break
          }
          Scalefree = sample_pa(n, power = p, m = EdgesEachStep, out.dist = NULL, out.seq = NULL,
                                out.pref = FALSE, zero.appeal = 1, directed = FALSE,
                                algorithm = "psumtree", start.graph = NULL)
          ##print(paste(" probablity ",p))
          # #print(paste(" vertices ",n))

          #calculation of properties from graph

          Degree = degree(Scalefree, normalized = FALSE)
          Closeness  = closeness.estimate(Scalefree, cutoff = -1)
          Estimate_betweenness  = estimate_betweenness(Scalefree, cutoff = -1)
          Mean_distance = mean_distance(Scalefree, directed = FALSE)
          Transitivity = transitivity(Scalefree, type = c("global"), vids = NULL, weights = NULL)
          Edges_get= gsize(Scalefree)
          # Post calculation of properties
          #print(paste(" Edges:",Edges_get))
          CentralityDegreeList[counter1]=list(Degree)

          CentralityClosenessList[counter1] = list(Closeness)

          CentralityBetweenessList[counter1] = list(Estimate_betweenness)

          CentralityDegreeMean[[counter1]] =  mean(Degree)
          CentralityClosenessMean[[counter1]] =  mean(Closeness)
          CentralityBetweenessMean[[counter1]] =  mean(Estimate_betweenness) #cutoff:The maximum path length to consider when calculating the betweenness. If zero or negative then there is no such limit.
          CentralityDegreeMin[[counter1]] =  min(Degree)
          CentralityClosenessMin[[counter1]] =  min(Closeness)
          CentralityBetweenessmin[[counter1]] =  min(Estimate_betweenness)
          CentralityDegreeMax[[counter1]] =  max(Degree)
          CentralityClosenessmax[[counter1]] =  max(Closeness)
          CentralityBetweenessmax[[counter1]] =  max(Estimate_betweenness)

          AvgGeodesicPath[[counter1]] = Mean_distance
          GlobalClusteringCoefficent[[counter1]] = Transitivity

          #Graph Info

          Vertices[[counter1]]<-n
          Probability[[counter1]]<-p
          NumberofEdges[[counter1]]<- Edges_get

          #Next graph counter increament

          counter1 = counter1+1;
          p= p+PstepSize;

        }

      }

    }
  }




  #return(CentralityDegreeList)

  return(list(data.frame(Vertices,NumberofEdges,Probability,GlobalClusteringCoefficent,CentralityDegreeMean,CentralityClosenessMean,CentralityBetweenessMean,AvgGeodesicPath),CentralityDegreeList,CentralityClosenessList));
  #return(data.frame(Vertices,Probability,GlobalClusteringCoefficent,CentralityDegreeMean,CentralityClosenessMean,CentralityBetweenessMean,AvgGeodesicPath, check.rows = FALSE));

}


SmallWorldGraphPropertiesEdges<- function(latticeDim,size1 = NULL,size2 = NULL,latticeNei,p1,p2,PstepSize,SizeVector = NULL,numberOfEdgesDelRandomly=NULL)
{
  #	Measure other properties of the graphs:
  #- Average geodesic/shortest path, global clustering coefficient, degree centrality(degree, closeness, betweeness)

  Probability<-1
  Vertices<-1
  NumberOfEdges<- 1

  CentralityDegreeMean<- 1
  CentralityClosenessMean<-1
  CentralityBetweenessMean<- 1
  CentralityDegreeMin<- 1
  CentralityClosenessMin<-1
  CentralityBetweenessmin<- 1
  CentralityDegreeMax<- 1
  CentralityClosenessmax<-1
  CentralityBetweenessmax<- 1

  CentralityDegreeList<-1
  names(CentralityDegreeList) = "Degree"
  CentralityClosenessList<-1
  names(CentralityClosenessList) = "Closeness"
  CentralityBetweenessList<-1
  names(CentralityBetweenessList) = "Estimate_betweenness"

  GlobalClusteringCoefficent <-1
  AvgGeodesicPath<- 1

  counter1 = 1
  if(!is.null(SizeVector))
  {
    for (n in SizeVector)
    {
      flag = 0;
      p=p1;

      repeat {

        if(p>p2)
        {
          break
        }

        smallWorldtest <- sample_smallworld(dim = latticeDim,size=n, nei=latticeNei, p=p)
        if(numberOfEdgesDelRandomly>0 && !is.null(numberOfEdgesDelRandomly))
        {
          for (delNum in 1:numberOfEdgesDelRandomly) {
            smallWorldtest = DelRandomEdge(smallWorldtest)
          }
        }

        #calculation of properties from graph
        ##print(paste(" probablity ",p))
        ##print(paste(" vertices ",n))
        Edges_get= gsize(smallWorldtest)
        #print(paste(" Edges:",Edges_get))

        Degree = degree(smallWorldtest, normalized = FALSE)
        Closeness  = closeness.estimate(smallWorldtest, cutoff = -1)
        Estimate_betweenness  = estimate_betweenness(smallWorldtest, cutoff = -1)
        Mean_distance = mean_distance(smallWorldtest, directed = FALSE)
        Transitivity = transitivity(smallWorldtest, type = c("global"), vids = NULL, weights = NULL)

        # Post calculation of properties

        CentralityDegreeList[counter1]=list(Degree)

        CentralityClosenessList[counter1] = list(Closeness)

        CentralityBetweenessList[counter1] = list(Estimate_betweenness)

        CentralityDegreeMean[[counter1]] =  mean(Degree)
        CentralityClosenessMean[[counter1]] =  mean(Closeness)
        CentralityBetweenessMean[[counter1]] =  mean(Estimate_betweenness) #cutoff:The maximum path length to consider when calculating the betweenness. If zero or negative then there is no such limit.
        CentralityDegreeMin[[counter1]] =  min(Degree)
        CentralityClosenessMin[[counter1]] =  min(Closeness)
        CentralityBetweenessmin[[counter1]] =  min(Estimate_betweenness)
        CentralityDegreeMax[[counter1]] =  max(Degree)
        CentralityClosenessmax[[counter1]] =  max(Closeness)
        CentralityBetweenessmax[[counter1]] =  max(Estimate_betweenness)

        AvgGeodesicPath[[counter1]] = Mean_distance
        GlobalClusteringCoefficent[[counter1]] = Transitivity

        #Graph Info

        Vertices[[counter1]]<-n
        Probability[[counter1]]<-p
        NumberOfEdges[[counter1]]<-Edges_get

        #Next graph counter increament

        counter1 = counter1+1;
        p= p+PstepSize;

      }

    }
  }
  else
  {
    for (n in size1:size2)
    {
      flag = 0;
      p=p1;

      repeat {

        if(p>p2)
        {
          break
        }

        smallWorldtest <- sample_smallworld(dim = latticeDim,size=n, nei=latticeNei, p=p)
        if(numberOfEdgesDelRandomly>0 && !is.null(numberOfEdgesDelRandomly))
        {
          for (delNum in 1:numberOfEdgesDelRandomly) {
            smallWorldtest = DelRandomEdge(smallWorldtest)
          }
        }

        #calculation of properties from graph
        ##print(paste(" probablity ",p))
        ##print(paste(" vertices ",n))
        Edges_get= gsize(smallWorldtest)
        ##print(paste(" Edges:",Edges_get))
        Degree = degree(smallWorldtest, normalized = FALSE)
        Closeness  = closeness.estimate(smallWorldtest, cutoff = -1)
        Estimate_betweenness  = estimate_betweenness(smallWorldtest, cutoff = -1)
        Mean_distance = mean_distance(smallWorldtest, directed = FALSE)
        Transitivity = transitivity(smallWorldtest, type = c("global"), vids = NULL, weights = NULL)

        # Post calculation of properties

        CentralityDegreeList[counter1]=list(Degree)

        CentralityClosenessList[counter1] = list(Closeness)

        CentralityBetweenessList[counter1] = list(Estimate_betweenness)

        CentralityDegreeMean[[counter1]] =  mean(Degree)
        CentralityClosenessMean[[counter1]] =  mean(Closeness)
        CentralityBetweenessMean[[counter1]] =  mean(Estimate_betweenness) #cutoff:The maximum path length to consider when calculating the betweenness. If zero or negative then there is no such limit.
        CentralityDegreeMin[[counter1]] =  min(Degree)
        CentralityClosenessMin[[counter1]] =  min(Closeness)
        CentralityBetweenessmin[[counter1]] =  min(Estimate_betweenness)
        CentralityDegreeMax[[counter1]] =  max(Degree)
        CentralityClosenessmax[[counter1]] =  max(Closeness)
        CentralityBetweenessmax[[counter1]] =  max(Estimate_betweenness)

        AvgGeodesicPath[[counter1]] = Mean_distance
        GlobalClusteringCoefficent[[counter1]] = Transitivity

        #Graph Info

        Vertices[[counter1]]<-n
        Probability[[counter1]]<-p
        NumberOfEdges[[counter1]]<-Edges_get

        #Next graph counter increament

        counter1 = counter1+1;
        p= p+PstepSize;

      }

    }
  }



  #return(CentralityDegreeList)

  return(list(data.frame(Vertices,Probability,NumberOfEdges,GlobalClusteringCoefficent,CentralityDegreeMean,CentralityClosenessMean,CentralityBetweenessMean,AvgGeodesicPath),CentralityDegreeList,CentralityClosenessList));
  #return(data.frame(Vertices,Probability,GlobalClusteringCoefficent,CentralityDegreeMean,CentralityClosenessMean,CentralityBetweenessMean,AvgGeodesicPath, check.rows = FALSE));



}



savePlot<- function(nameX,nameY,plotType,ParamFileNameDir)
{
  ParamFileName = ParamFileNameDir


  XnameString = gsub("[[:space:]]", "", nameX)
  YnameString = gsub("[[:space:]]", "", nameY)
  XnameString = str_replace_all(XnameString, "`", "")
  YnameString = str_replace_all(YnameString, "`", "")

  fileName<- paste(XnameString,YnameString,plotType,".pdf",sep ="")
  ggsave(paste(ParamFileName,"/",ParamFileName,fileName,sep = ""), width = 20, height = 20, units = "cm")


}



Plot2dListOfRandomGraphPropertiesMean<- function(givenObject,x,y,xlabel,ylabel)
{

  #This functions plots the average value from the multiple graph generated by the same parameteres
  totalNumberofSamples = length(givenObject)

  if(is.null(dim(givenObject[[1]]))==FALSE){
    Samples = givenObject
    xCounter<-0
    yCounter<-0
    counter<-0

    for (Sample in Samples) {

      xCounter = Sample[[x]]+xCounter
      yCounter = Sample[[y]]+yCounter
      counter=counter+1

    }

    # colNumbers = ncol(Samples[[1]])

    xMean = xCounter/counter
    yMean = yCounter/counter



    # plot3d(xMean, yMean,xlab = xlabel,ylab = ylabel)


    #surf3D(xMean, yMean, zMean, phi = 45, theta = 45,xlab = xlabel,ylab = ylabel,zlab = zlabel)




    df = data.frame(xMean,yMean)

    names(df)[names(df) == 'xMean'] <- xlabel
    names(df)[names(df) == 'yMean'] <- ylabel


    cat("test")


    return(df)

  }
  else if(is.null(dim(givenObject))==TRUE){

    xCounter<-0
    yCounter<-0
    zCounter<-0
    counter<-0

    for (i in 1:totalNumberofSamples) {
      Sample=as.data.frame(givenObject[[i]][[1]])

      xCounter = Sample[[x]]+xCounter
      yCounter = Sample[[y]]+yCounter

      counter=counter+1

    }

    # colNumbers = ncol(Samples[[1]])

    xMean = xCounter/counter
    yMean = yCounter/counter



    # plot(xMean, yMean,xlab = xlabel,ylab = ylabel)


    #surf3D(xMean, yMean, zMean, phi = 45, theta = 45,xlab = xlabel,ylab = ylabel,zlab = zlabel)




    df = data.frame(xMean,yMean)

    names(df)[names(df) == 'xMean'] <- xlabel
    names(df)[names(df) == 'yMean'] <- ylabel


    # cat("test")


    return(df)


  }

}


