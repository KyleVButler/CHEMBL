
probe_list <- read_csv("PROBELIST.csv")
#----------------------------------------------
#get first and second degree connections to target, find all probes, and all nodes that connect target to probe
network_table <- read_csv("network_table.csv")
network_table <- filter(network_table, confidenceScore >= 0.25)
target <- "P01116"
probe_genes <- as.character(unique(probe_list$UNIPROTKB))
#target_connections <- network_table[apply(network_table, 1, FUN = 
#                                 function(x, y = target) {any(grepl(y, x, ignore.case = TRUE))}), ]
target_connections <- network_table[network_table$A == target | network_table$B == target, ]
target_primary <- unique(c(target_connections$A, target_connections$B))
target_connections <- network_table[network_table$A %in% target_primary | network_table$B %in% target_primary, ]
relevant_probes <- probe_genes[probe_genes %in% unique(c(target_connections$A, target_connections$B))]
probe_primary <- network_table[network_table$A %in% relevant_probes | network_table$B %in% relevant_probes, ]
probe_primary <- unique(c(probe_primary$A, probe_primary$B))
middle_nodes <- probe_primary[probe_primary %in% target_primary]
#build attribute list
name_list <- tibble(node = c(target_connections$A, target_connections$B), alias = c(target_connections$aliasA, target_connections$aliasB))
name_list$char <- nchar(name_list$alias)
name_list <- name_list %>% group_by(node) %>% arrange(char) %>% dplyr::slice(1)
name_list$char[is.na(name_list$char)] <- 0
name_list$alias[name_list$char > 8] <- name_list$node[name_list$char > 8]
name_list$alias[name_list$char == 0] <- name_list$node[name_list$char == 0]
relevant_nodes <- tibble(node = unique(c(middle_nodes, target, relevant_probes)), node_type = "middle")
relevant_nodes <- left_join(relevant_nodes, name_list)
relevant_nodes[relevant_nodes$node == target, 2] <- "target"
relevant_nodes[relevant_nodes$node %in% relevant_probes, 2] <- "probe"


target_connections <- network_table[network_table$A %in% relevant_nodes$node & network_table$B %in% relevant_nodes$node, ]
target_connections$interaction_type <- "direct"
target_connections$interaction_type[grepl("reaction", target_connections$type, ignore.case = TRUE)] <- "reaction"
g <- graph.data.frame(as.matrix(target_connections[, c(1, 2, 8)]), vertices = relevant_nodes, directed = FALSE)
E(g)$color <- ifelse(E(g)$interaction_type == "reaction", "orange", "grey")
V(g)[V(g)$node_type == "probe"]$color <- "red"
V(g)[V(g)$node_type == "middle"]$color <- "yellow"
V(g)[V(g)$node_type == "target"]$color <- "lightblue"
V(g)$name <- V(g)$alias
plot(g, layout = layout.fruchterman.reingold, vertex.size = 4, vertex.label.cex = 0.9, vertex.label.dist=0.35, edge.curved = FALSE)

