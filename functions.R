### DrawDAG

DrawDAG = function(adjacency, layers, prop = NULL, edge.col = "#000000", edge.alpha = 0.5,
                   node.col = "#ffffff", node.size = 5, direct = FALSE){
  ### Graph
  mygraph = Graph(adjacency = adjacency, mode = "directed")
  
  has.path = (apply(adjacency,1,sum)+apply(adjacency,2,sum))!=0
  
  ### Edge ID
  start_id = ends(mygraph, es = E(mygraph))[,1] # start point
  end_id = ends(mygraph, es = E(mygraph))[,2] # end point
  
  if(is.null(prop)){
    myedge_colours = edge.col[1]
  } else {
    ### Edge colour
    myedge_colours = ifelse(posprop[adjacency==1]>0.5, edge.col[1],
                            ifelse(posprop[adjacency==1]<0.5, edge.col[2],edge.col[3]))
    names(myedge_colours) = paste0(start_id,end_id)
  }
  myedge_colours = alpha(myedge_colours, edge.alpha)

  ### Vertex ID
  vertex_id = colnames(adjacency)[has.path]
  outcome_id = vertex_id[vcount(mygraph)]
  
  ### Vertex colour
  mynode_colours = node.col[has.path]
  tmp = adjacency[has.path, has.path]
  
  V(mygraph)$color = mynode_colours
  V(mygraph)$frame.color = ifelse(tmp[,ncol(tmp)]==1, myedge_colours[paste0(vertex_id,outcome_id)],NA)
  
  # #### Test concordance
  # node.sign = V(mygraph)$frame.color
  # names(node.sign) = vertex_id
  # 
  # sign = ifelse(myedge_colours != node.sign[end_id], 1, -1) # +/- and +/- are concordant associations for protective factors
  # names(sign) = start_id
  # 
  # risk = colnames(dag)[c(1,4)]
  # sign[start_id %in% risk] = sign[start_id %in% risk]*-1
  # myedge_colours = ifelse(end_id != outcome_id & sign==1, "#000000", myedge_colours)
  
  ### Vertex size
  V(mygraph)$size = node.size[has.path]
  
  ### Vertex shape
  add_shape("nil")
  V(mygraph)$shape = rep(c("square","circle","nil"),pk)[has.path]
  
  ### Vertex label
  V(mygraph)$label = c(colnames(expo),mylabels[colnames(X2)], "BMI")[has.path]
  
  if(direct){
    myedge_colours = ifelse(names(myedge_colours) %in% paste0(vertex_id,outcome_id), myedge_colours, NA)
    V(mygraph)$size = ifelse(tmp[,ncol(tmp)]==1, V(mygraph)$size, 0)
    V(mygraph)$label = ifelse(tmp[,ncol(tmp)]==1, V(mygraph)$label, "")
    V(mygraph)$label[vcount(mygraph)] = "BMI"
  }
  # Edge
  E(mygraph)$color = myedge_colours
  E(mygraph)$width = 2
  
  # E(mygraph)$width = 2*selprop[adjacency==1]
  
  # #### Transparency (effect magnitude)
  # weights = abs(betas[adjacency==1])
  # E(mygraph)$color = alpha(myedge_colours, 0.1 + log(weights/max(weights) + 1)) # Log-transform the normalized weights
  
  # Arrow
  E(mygraph)$arrow.size = 1
  E(mygraph)$arrow.width = 1
  E(mygraph)$arrow.mode = "->"
  
  direct_link = incident(mygraph, vcount(mygraph), mode = "in")
  E(mygraph)$arrow.mode[direct_link] = "-"
  
  ### Layout
  layers = layers[has.path]
  layout = layout_with_sugiyama(mygraph, layers = layers)
  
  set.seed(1)
  plot(mygraph, layout = layout, asp = asp, vertex.label.cex = 0.7)
  legend("bottomleft",
         legend = c(paste0(c("Positive","Negative")," associations")),
         bty = "n", lty = rep(1,3),
         col = c("#dc251f","#0019a8"), cex = 0.6)
  # legend("bottomleft",
  #        legend = c(paste0(c("Risk","Protective")," biomarkers"), paste0(c("Positive","Negative")," associations")),
  #        pt.cex = c(rep(1.2,2),rep(NA,2)),
  #        pch = c(rep(1,2),rep(NA,2)),
  #        bty = "n", lty = c(rep(NA,2),rep(1,3)),
  #        col = c("#dc251f","#0019a8"), cex = 0.6)
}
