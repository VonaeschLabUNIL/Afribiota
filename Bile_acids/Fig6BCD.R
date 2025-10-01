
############### Workspace preparation #####################################################################

# load libraries
pacman::p_load(readr, ggplot2, dplyr, ggpubr, tidyr, forcats)

input_path = "/Users/maier/Desktop/Code for Paper/Data/"
output_path = "/Users/maier/Desktop/Code for Paper/Output/"

bile_acids <- read_tsv(paste0(input_path,"bile_acids_monocultures.tsv" ))


############### Plot data  #####################################################################

order_species <- bile_acids %>%  group_by(name) %>% summarize(sum = sum(hit)) %>% arrange(desc(sum))
order_species <- as.vector(order_species$name)
bile_acids <- bile_acids %>% mutate(name=factor(name, levels = order_species))

bile_acids <- bile_acids %>% mutate(category_drug_2=factor(category_drug_2, levels=c('7-keto Lithocholic acid', 'Deoxycholic acid', 'Lithocholic acid', 'Lithoderivatives', 'Hyodeoxycholic acid', 'Dehydrolithocholic acid','NA')))
bile_acids <- bile_acids %>% mutate(category_drug=factor(category_drug, levels=c('primary', 'conjugated primary', 'secondary', 'conjugated secondary', 'other')))


myplot <-ggplot(bile_acids, aes(x = name, y = drug_name_3)) +
  geom_tile(aes(fill = AUC)) +
  theme_bw() +
  scale_fill_gradient2(low = "#fc8d59",mid = "#ffffbf",high = "#91cf60",midpoint = 1,space = "Lab",na.value = "grey50")+
  theme(strip.background = element_blank(),axis.title.x=element_blank(), axis.title.y=element_blank(),axis.text.y= element_text(size=6), axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))+
  facet_grid(category_drug~category, scales='free',space = 'free')

ggsave(filename=paste0(output_path,"Fig6D.pdf"), plot = myplot, device = cairo_pdf, width = 297, height = 210, units = "mm")

hit_per_category <- bile_acids %>%  group_by(category, category_drug) %>% summarize(sum = sum(hit)) %>% arrange(desc(sum))
hit_per_category <- bile_acids %>%  group_by(category, category_drug) %>% count(hit) 
hit_per_category <- hit_per_category %>%  spread(hit, n)
names(hit_per_category)[names(hit_per_category)=="TRUE"] <- "nr_hits"
names(hit_per_category)[names(hit_per_category)=="FALSE"] <- "nr_none_hits"

hit_per_category <- hit_per_category %>% mutate(all = nr_hits +  nr_none_hits)
hit_per_category <- hit_per_category %>% mutate(percent_hits = nr_hits/all*100)
hit_per_category <- hit_per_category %>% mutate(label = paste0(round(percent_hits,1), '% (', nr_hits, '/', all, ')'))

myplot2 <- ggplot(hit_per_category, aes(x=fct_rev(category_drug), y=percent_hits, fill=fct_rev(category))) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values=c("#fc8d59","#91cf60"))+
  geom_text( aes(x=fct_rev(category_drug), y=45, label=label, group=fct_rev(category)), position = position_dodge(width = .9),size = 3) +
  theme_minimal()+
  theme(axis.text.x = element_text(size = 14),  axis.text.y = element_text(size = 14), legend.title = element_blank(),  axis.title.x = element_blank(),  axis.title.y = element_blank())+
  coord_flip() +
  ggtitle('Inhibitory effect of bile acids')

ggsave(filename=paste0(output_path,"Fig6B.pdf"), plot = myplot2, device = cairo_pdf, width = 297, height = 210, units = "mm")


hit_per_category2 <- bile_acids %>%  group_by(category, category_drug_2) %>% count(hit) 
hit_per_category2 <- hit_per_category2 %>%  spread(hit, n)
names(hit_per_category2)[names(hit_per_category2)=="TRUE"] <- "nr_hits"
names(hit_per_category2)[names(hit_per_category2)=="FALSE"] <- "nr_none_hits"

hit_per_category2 <- hit_per_category2 %>% mutate(all = nr_hits +  nr_none_hits)
hit_per_category2 <- hit_per_category2 %>% mutate(percent_hits = nr_hits/all*100)
hit_per_category2 <- hit_per_category2 %>% mutate(label = paste0(round(percent_hits,1), '% (', nr_hits, '/', all, ')'))

myplot5 <-  ggplot(filter(bile_acids, category_drug_2 != 'NA'), aes(x=fct_rev(category_drug_2), y=AUC, fill=fct_rev(Lachnospiraceae))) + 
  geom_boxplot(position=position_dodge(1)) +
  theme_minimal()+
  coord_flip() +
  scale_fill_manual(values=c("#91cf60","#fc8d59"))+
  theme(axis.text.x = element_text(size = 14),  axis.text.y = element_text(size = 14), legend.title = element_blank(),  axis.title.x = element_blank(),  axis.title.y = element_blank())

ggsave(filename=paste0(output_path,"Fig6C.pdf"), plot = myplot5, device = cairo_pdf, width = 297, height = 210, units = "mm")

