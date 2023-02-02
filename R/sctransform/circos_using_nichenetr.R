
d2 <- DotPlot(seuratObj, features = best_upstream_ligands, split.by = "type") + RotatedAxis()

png(file="top_ligands_dotplot.png",
    width=700, height=1000)
d1/d2
dev.off()

avg_expression_ligands = AverageExpression(seuratObj, features =best_upstream_ligands)

sender_ligand_assignment = avg_expression_ligands$RNA %>% apply(1, function(ligand_expression){
  ligand_expression > (ligand_expression %>% mean() + ligand_expression %>% sd())
}) %>% t()
sender_ligand_assignment = sender_ligand_assignment %>% apply(2, function(x){x[x == TRUE]}) %>% purrr::keep(function(x){length(x) > 0})
names(sender_ligand_assignment)

all_assigned_ligands = sender_ligand_assignment %>% lapply(function(x){names(x)}) %>% unlist()
unique_ligands = all_assigned_ligands %>% table() %>% .[. == 1] %>% names()
general_ligands = best_upstream_ligands %>% setdiff(unique_ligands)

Fibroblasts3_specific_ligands = sender_ligand_assignment$'Fibroblasts 3' %>% names() %>% setdiff(general_ligands)
Fibroblasts4_specific_ligands = sender_ligand_assignment$'Fibroblasts 4' %>% names() %>% setdiff(general_ligands)
Fibroblasts5_specific_ligands = sender_ligand_assignment$'Fibroblasts 5' %>% names() %>% setdiff(general_ligands)
Fibroblasts6_specific_ligands = sender_ligand_assignment$'Fibroblasts 6' %>% names() %>% setdiff(general_ligands)
Fibroblasts7_specific_ligands = sender_ligand_assignment$'Fibroblasts 7' %>% names() %>% setdiff(general_ligands)
Macrophages1_specific_ligands = sender_ligand_assignment$'Macrophages 1' %>% names() %>% setdiff(general_ligands)
Macrophages2_specific_ligands = sender_ligand_assignment$'Macrophages 2' %>% names() %>% setdiff(general_ligands)
Sebaceousglandcells1_specific_ligands = sender_ligand_assignment$'Sebaceous gland cells 1' %>% names() %>% setdiff(general_ligands)
Sebaceousglandcells2_specific_ligands = sender_ligand_assignment$'Sebaceous gland cells 2' %>% names() %>% setdiff(general_ligands)
Endothelial2_specific_ligands = sender_ligand_assignment$'Endothelial 2' %>% names() %>% setdiff(general_ligands)
SkMuscle1_specific_ligands = sender_ligand_assignment$'Sk Muscle 1' %>% names() %>% setdiff(general_ligands)
Adipocyte_specific_ligands = sender_ligand_assignment$'Adipocyte' %>% names() %>% setdiff(general_ligands)
Immune_specific_ligands = sender_ligand_assignment$'Immune' %>% names() %>% setdiff(general_ligands)


ligand_type_indication_df = tibble(
  ligand_type = c(rep("Fib 3", times = Fibroblasts3_specific_ligands %>% length()),
                  rep("Fib 4", times = Fibroblasts4_specific_ligands %>% length()),
                  rep("Fib 5", times = Fibroblasts5_specific_ligands %>% length()),
                  rep("Fib 6", times = Fibroblasts6_specific_ligands %>% length()),
                  rep("Fib 7", times = Fibroblasts7_specific_ligands %>% length()),
                  rep("Macro 1-specific", times = Macrophages1_specific_ligands %>% length()),
                  rep("Macro 2-specific", times = Macrophages2_specific_ligands %>% length()),
                  rep("SBG 1-specific", times = Sebaceousglandcells1_specific_ligands %>% length()),
                  rep("SBG 2-specific", times = Sebaceousglandcells2_specific_ligands %>% length()),
                  rep("Endo 2-specific", times = Endothelial2_specific_ligands %>% length()),
                  rep("SK Mus 1-specific", times = SkMuscle1_specific_ligands %>% length()),
                  rep("Adipo-specific", times = Adipocyte_specific_ligands %>% length()),
                  rep("Imm-specific", times = Immune_specific_ligands %>% length()),
                  rep("General", times = general_ligands %>% length())),
  ligand = c(Fibroblasts3_specific_ligands,Fibroblasts4_specific_ligands,Fibroblasts5_specific_ligands,Fibroblasts6_specific_ligands,Fibroblasts7_specific_ligands,
             Macrophages1_specific_ligands,Macrophages2_specific_ligands,Sebaceousglandcells1_specific_ligands,Sebaceousglandcells2_specific_ligands,
             Endothelial2_specific_ligands,SkMuscle1_specific_ligands,Adipocyte_specific_ligands,Immune_specific_ligands,general_ligands))
active_ligand_target_links_df = active_ligand_target_links_df %>% mutate(target_type = "Injured-DE") %>% inner_join(ligand_type_indication_df)
cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(0.40)
active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)

ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())

circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)

grid_col_ligand =c("General" = "lawngreen",
                   "Fib 4" = "royalblue",
                   "Macro 1-specific" = "darkgreen",
                   "Adipo-specific" = "violet",
                   "Fib 3" = "steelblue2",
                   "Fib 6" = "darkblue",
                   "Fib 5" ="blueviolet",
                   "Fib 7" ="blue3",
                   "SK Mus 1-specific" = "yellow"
                   
                   )
grid_col_target =c(
  "Injured-DE" = "tomato")

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target = tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
links_circle = circos_links %>% select(ligand,target, weight)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(target,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$target)

grid_col =c(grid_ligand_color,grid_target_color)

# give the option that links in the circos plot will be transparant ~ ligand-target potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 
target_order = circos_links$target %>% unique()
ligand_order = c(Fibroblasts3_specific_ligands,Fibroblasts4_specific_ligands,Fibroblasts5_specific_ligands,Fibroblasts6_specific_ligands,Fibroblasts7_specific_ligands,
                 Macrophages1_specific_ligands,Adipocyte_specific_ligands,SkMuscle1_specific_ligands, general_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,target_order)

width_same_cell_same_ligand_type = 1
width_different_cell = 10
width_ligand_target = 15
width_same_cell_same_target_type = 1

gaps = c(
  # width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Fib 3") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Fib 4") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Fib 5") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Fib 6") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Fib 7") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Macro 1-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Adipo-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "SK Mus 1-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "General") %>% distinct(ligand) %>% nrow() -1)),
  width_ligand_target,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_type == "Injured-DE") %>% distinct(target) %>% nrow() -1)),
  width_ligand_target
)
library(circlize)
#gaps <- c(-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1)
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = 0, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA)
circos.clear()

svg("ligand_receptor_circos.svg", width = 15, height = 15)
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #

