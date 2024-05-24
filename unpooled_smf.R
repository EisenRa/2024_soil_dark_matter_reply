mf <- read_delim("data/microbial_fractions.csv.gz")
md <- read_delim("data/metadata.tsv") %>%
  clean_names()
enviro_data <- readxl::read_excel("data/Global_Data.xlsx", sheet = 7) %>%
  clean_names()
avg_genome_size <- readxl::read_excel("data/Global_Data.xlsx", sheet = 1) %>% 
  select(sample_name, `Average Genome size (bp)`) %>%
  clean_names()

df <- mf %>%
  inner_join(md, by = join_by("sample" == "run"), suffix=c("", ".y")) %>%
  select(-ends_with(".y")) %>% 
  inner_join(enviro_data, by = join_by(sample_name)) %>% 
  inner_join(avg_genome_size, by = join_by(sample_name)) %>%
  clean_names()

df %>% 
  summarise(n = n(), .by = sample_name) %>% 
  arrange(desc(n))


#Immediately, it's a bit confusing, as the bioproject PRJEB18701 seems to contain more (288) than the 128 samples specified by Piton et al. 2023:
##"We analysed a global dataset of 128 metagenomes each from unique soil samples distributed across continents and latitude (Extended Data Fig. 8)."

#It looks like there were multiple sequencinig runs of the sample samples, so will merge them by biosample.
#The SMF (SingleM Microbial Fraction) values for accessions with the same sample_name are all within <5%, so I'll assume this is the case, and take the average:

df_filt <- df %>%
  summarise(smf = mean(read_fraction),
            AGS_smf = mean(average_bacterial_archaeal_genome_size),
            .by = sample_name) %>%
  left_join(df, by = join_by(sample_name)) %>%
  distinct(sample_name, .keep_all = TRUE) %>%
  type.convert(., as.is =TRUE)  

#We're left now with 126 samples (close to the original 128).


df_filt %>%
  ggplot(aes(x = environment_biome, y = (100 - smf), 
             colour = environment_biome)) +
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.3, alpha = 0.7) +
  geom_hline(yintercept = (100 - mean(df_filt$smf)),
             linetype = 2) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  ylab("Non-bacterial fraction (%)")


df_filt %>%
  ggplot(aes(x = soilp_h, y = AGS_smf / 1000000)) +
  geom_smooth(method = "lm") +
  geom_point(aes(colour = environment_biome)) +
  theme_classic() +
  labs(x = "Soil pH", y = "SingleM Average Genome Size (Mbp)")

cor(df_filt$AGS_smf, df_filt$soilp_h)