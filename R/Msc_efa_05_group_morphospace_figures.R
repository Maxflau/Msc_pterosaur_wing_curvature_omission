# Save the remaining group morphospace figures for the last categories.
save_plot(make_group_morph("Environment",
                           c("Archipelago"="#A300E9","Coastal environment"="#4F4F4F",
                             "Desert"="#B51963","Floodplain"="#FFC107",
                             "Fluvial environment"="#FF5722","Forested wetland environment"="#FF9800",
                             "Marine environment"="#FAAF90","Semi-arid floodplain"="#4CAF50"),
                           "Morphospace Occupation by Environment"),"CLADE_11_morphospace_environment_convex_hulls")

save_plot(make_group_morph("Diet.1",
                           c("Carnivore"="#E74C3C","Durophageous"="#4CAF50","Filter Feeder"="#7B1FA2",
                             "Generalist"="#144B00","Herbivore"="#80DEEA","Insectivore"="#FFA726",
                             "Piscivore"="#F5F5B8"),
                           "Morphospace Occupation by Primary Diet"),
          "CLADE_12_morphospace_diet_convex_hulls")

save_plot(make_group_morph("Diet_combined",
                           c("Carnivore+Durophageous"  ="#9575CD","Carnivore+Generalist"    ="#B0BEC5","Carnivore+Insectivore"   ="#42A5F5",
                             "Carnivore+Piscivore"     ="#7B1FA2","Durophageous+Piscivore"  ="#FF8A65","Filter Feeder+Piscivore" ="#FF7043",
                             "Generalist+Piscivore"    ="#EC407A","Herbivore+Durophageous"  ="#F48FB1","Insectivore+Carnivore"   ="#26A69A",
                             "Insectivore+Durophageous"="#CE93D8","Insectivore+Piscivore"   ="#EF5350","Piscivore+Carnivore"     ="#FFA726",
                             "Piscivore+Durophageous"  ="#66BB6A","Piscivore+Generalist"    ="#FFEE58","Piscivore+Insectivore"   ="#FFC107",
                             "Piscivore+Piscivore"     ="#FDD835"),
                           "Morphospace by Diet Combination",
                           "Convex hulls by combined diet") +
            # Legend style matching Image 3: large circles, bold title "Diet (combined)"
            guides(fill   = guide_none(),colour = guide_legend(title="Diet (combined)", ncol=1,override.aes=list(shape=21, size=4.5, fill=NA, stroke=1.0))) +
            theme(legend.text  = element_text(size=8.5),
                  legend.title = element_text(size=11, face="bold"),
                  plot.margin  = unit(c(8, 190, 8, 8), "pt")),
          "CLADE_13_morphospace_dietcombo_convex_hulls", w=14
)
save_plot(
  make_group_morph("Depositional.settings.paleoenvironment",
                   c("Aeolian"="#482072","Alluvial plain"="#4CAF50","Coastal/Shallow marine"="#D55E00",
                     "Fluviodeltaic"="#144B00","Lacustrine: large lake"="#FFC107","Lacustrine: small lake"="#F0CD5D",
                     "Lagoonal deposit"="#A300E9","playa"="#6A9DA3"),"Morphospace Occupation by Depositional Setting",
                   "Taphonomic context and wing functional diversity"),
  "CLADE_14_morphospace_depositional_convex_hulls")

cat("\nAll morphospace figures saved.\n")
