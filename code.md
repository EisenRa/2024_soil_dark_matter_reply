# Re-analysis of soil replies with SMF values
Raphael Eisenhofer

### Import and clean data

``` r
library(tidyverse)
library(janitor)

pooled_smf <- read_delim("data/smf.csv")
md <- read_delim("data/metadata.tsv") %>%
  clean_names()
enviro_data <- readxl::read_excel("data/Global_Data.xlsx", 
                                  sheet = 7) %>%
  clean_names()
avg_genome_size <- readxl::read_excel("data/Global_Data.xlsx", 
                                      sheet = 1) %>% 
  select(sample_name, `Average Genome size (bp)`) %>%
  clean_names()

# join tables
df <- pooled_smf %>%
  inner_join(md, by = join_by("sample" == "sample_alias"), 
             suffix=c("", ".y")) %>%
  select(-ends_with(".y")) %>% 
  inner_join(enviro_data, by = join_by(sample_name)) %>% 
  inner_join(avg_genome_size, by = join_by(sample_name)) %>%
  distinct(sample, .keep_all = T) %>%
  type.convert(., as.is =TRUE) %>%
  clean_names()
```

We’re left now with 126 samples (close to the original 128).

### How much non-bacterial DNA are in these samples?

Piton et al. asserted that:

> “However, given that eukaryotic sequences represent \<2% of annotated
> sequences, we might also expect only a small fraction of eukaryotic
> sequences in the unannotated base pairs.”
>
> “Indeed, assuming an extreme range of 4–9% eukaryotic base pairs..”

Let’s now use our SMF (**S**ingleM **M**icrobial **F**raction) estimates
to check how much non-bacterial DNA there is in these samples:

``` r
df %>%
  ggplot(aes(x = environment_biome, y = (100 - read_fraction), 
             colour = environment_biome)) +
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.3, alpha = 0.7) +
  geom_hline(yintercept = (100 - mean(df$read_fraction)),
             linetype = 2) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
  ylab("Non-bacterial fraction (%)")
```

![](code_files/figure-commonmark/unnamed-chunk-2-1.png)

OK, so the assertion that 4-9% of the base pairs are eukaryotic is a
severe underestimate! 10.1- to 4.5 -fold understimations to be precise.

It’s also interesting that we do see variation in non-bacterial read
fractions between soil environmental biomes.

### How accurate are Osmund and Piton et al.’s AGS estimates?

Piton et al.’s original AGS (**A**verage **G**enome **S**ize) estimates
are based on the tool MicrobeCensus. Here’s the summary from the
[original 2015
paper](https://link.springer.com/article/10.1186/s13059-015-0611-7):

> “Specifically, the AGS of a community will be inversely proportional
> to the relative abundance, R, of an essential single-copy gene in that
> community: AGS ∝ R -1. In other words, these essential genes will be
> sequenced at a higher rate in a community with a small AGS relative to
> a community with a large AGS; this is simply because these genes make
> up a larger fraction of the total genomic DNA in the community with
> smaller genomes.”

So essentially, a major assumption here is that most (if not all) of the
metagenome is bacterial. If you had substantial euk DNA, you would
severely overestimate the AGS of the metagenome — as the single-copy
core genes will be artificially lower in relative abundance -\> higher
AGS estimate.

Osmund et al. were correct in pointing out this major assumption. In
their reply, they calculated AGS using only reads that aligned to
bacterial genomes, and reported a much lower AGS: 3 Mbp, vs Piton et
al.’s 6.8 Mbp.

SMF can help here, as it:

1.  Accounts for the fact that different metagenomes can have different
    proportions of non-bacterial DNA
2.  Is not reliant on alignment to reference genomes

Let’s see where SMF’s AGS estimates lie between Piton et al. and Osmund
et al.s’:

``` r
SMF_AGS <- as.numeric(scales::number(mean(df$average_bacterial_archaeal_genome_size) / 1000000, accuracy = 0.1))

AGS_overestimate_fold <- as.numeric(scales::number(6.8 / SMF_AGS, accuracy = 0.01))
AGS_overestimate_pc <- (1 - as.numeric(scales::number(SMF_AGS / 6.8, accuracy = 0.001))) * 100
```

So SMF’s AGS estimate 4.7 Mbp is between Osmund et al.’s 3 Mbp and Piton
et al.’s 6.8 Mbp. With Piton et al.’s 6.8 Mbp being a 1.45 -fold
overestimate (30.9 %).

### Does Piton et al.’s original relationship between AGS and soil pH still supported?

Now that we have more robust AGS estimates thanks to SMF, let’s see if
this relationship is still supported. This was Osmund et al.’s main
point:

> “Therefore, we suggest that the strong association of very large AGS
> with acidic soils in the full metagenomes is likely an artefact of
> ecosystems with acidic soils having larger proportions of
> non-bacterial DNA.”

``` r
df %>%
  ggplot(aes(x = soilp_h, y = average_bacterial_archaeal_genome_size / 1000000)) +
  geom_smooth(method = "lm") +
  geom_point(aes(colour = environment_biome)) +
  theme_classic() +
  labs(x = "Soil pH", y = "SingleM Average Genome Size (Mbp)")
```

    `geom_smooth()` using formula = 'y ~ x'

![](code_files/figure-commonmark/unnamed-chunk-4-1.png)

``` r
correlation <- cor(df$average_bacterial_archaeal_genome_size, df$soilp_h)
```

Despite their inaccurate AGS estimates, we still see a strong negative
correlation as observed by Piton et al. Correlation = -0.7633485.
