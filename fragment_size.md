# mpcrseq exploratory

Import data and join amplicon size data with genotyping data.


```r
library(tidyverse)
library(lme4)

vcf	               <- read_csv("vcf_129snp_long.csv")
```

```
## Parsed with column specification:
## cols(
##   CHROM = col_double(),
##   POS = col_double(),
##   REF = col_character(),
##   ALT = col_character(),
##   SAMPLE_NAME_IN_RUN = col_character(),
##   GT = col_character(),
##   PL = col_character(),
##   DP = col_double(),
##   AD = col_character(),
##   GQ = col_double()
## )
```

```r
s1s2_product_sizes <- read_tsv("set12_product_sizes.tsv")
```

```
## Parsed with column specification:
## cols(
##   snp_id = col_character(),
##   template_size = col_double(),
##   pcr1_product = col_double(),
##   pcr2_product = col_double()
## )
```

```r
s1s2$snp_id <- gsub('-500','', s1s2$snp_id)

vcf$snp_id <- paste(vcf$CHROM, vcf$POS, sep = ":")

vcf <- left_join(s1s2, vcf)
```

```
## Joining, by = "snp_id"
```

```r
vcf <- vcf %>% filter(!is.na(CHROM))

vcf$sample <- vcf$SAMPLE_NAME_IN_RUN

metadata <- read_csv('big_bam_list.csv')
```

```
## Parsed with column specification:
## cols(
##   sample = col_character(),
##   run_id = col_character(),
##   sample_id = col_character(),
##   species = col_character(),
##   type = col_character(),
##   dna = col_double(),
##   Year = col_character(),
##   region = col_character()
## )
```

```r
vcf <- left_join(vcf, metadata)
```

```
## Joining, by = "sample"
```


Generate some variables, and group data by snp_id.


```r
# Separate allele depth & genotype likelihoods

vcf <- vcf %>% separate(AD, c("AD1", "AD2"), ",", convert = TRUE, remove = FALSE) %>% separate(PL, c("pl00", "pl01", "pl11"), convert = TRUE, remove = FALSE)
```

```
## Warning: Expected 2 pieces. Missing pieces filled with `NA` in 600 rows [3901, 3902, 3903, 3904, 3905, 3906, 3907, 3908, 3909, 3910, 3911, 3912, 3913, 3914, 3915, 3916, 3917, 3918, 3919, 3920, ...].
```

```
## Warning: Expected 3 pieces. Missing pieces filled with `NA` in 600 rows [3901, 3902, 3903, 3904, 3905, 3906, 3907, 3908, 3909, 3910, 3911, 3912, 3913, 3914, 3915, 3916, 3917, 3918, 3919, 3920, ...].
```

```r
# Calc. likelihood ratio
likelihood_ratio <- function(x) {
  y <- 10^(-x/10)
  max(y) / (sum(y) - max(y))
}

vcf <- vcf %>% rowwise() %>% mutate(lr = likelihood_ratio(c(pl00,pl01,pl11)))

# Group by snp_id, calc snp-level variables.

snp_data <- vcf %>% 
  group_by(snp_id, CHROM, REF, ALT, template_size, pcr1_product, pcr2_product) %>% 
  summarize(n_called = sum(called), n = n(), call_rate = n_called/n, AD1 = sum(AD1),
	          AD2 = sum(AD2), totalDP = sum(DP), meanDP = mean(DP), medianDP = median(DP),
						paa = sum(gt == 0) / n(), pAa = sum(gt == 1) / n(), pAA = sum(gt == 2) / n())
```

```
## Warning: Grouping rowwise data frame strips rowwise nature
```

```
## Error: object 'called' not found
```

```r
# Minor allele frequency

getmaf <- function(ad1, ad2) {
  if(ad1 == 0 || is.na(ad1) || ad2 == 0 || is.na(ad2)) {return(0)}

  af1 <- ad1 / (ad1 + ad2)
  if(af1 > 0.5){
    return(1 - af1)
  } else {
    return(af1)
  }
}

snp_data <- snp_data %>% rowwise() %>% mutate(maf = getmaf(AD1, AD2))

# Allele frequency
getaf <- function(ad1, ad2) {

  if(ad2 == 0 || is.na(ad2)) {return(0)}
	if(ad1 == 0 || is.na(ad1)) {return(1)}

  af1 <- ad1 / (ad1 + ad2)
}

snp_data <- snp_data %>% rowwise() %>% mutate(af = getaf(AD1, AD2))



sample_data <- vcf %>% group_by(sample, type, dna, Year, region, run_id) %>% summarise(n_called = sum(called), n = n(), call_rate = n_called/n, AD1 = sum(AD1),
	          AD2 = sum(AD2), totalDP = sum(DP), meanDP = mean(DP), medianDP = median(DP))
```

```
## Warning: Grouping rowwise data frame strips rowwise nature
```

```
## Error: object 'called' not found
```

```r
sample_data <- sample_data %>% group_by(run_id) %>% mutate(totalRunDP = sum(totalDP)) %>% mutate(scaledDP = totalDP / totalRunDP)
```

SNP Depth


```r
vcf$SNPID <- factor(vcf$snp_id, levels = snp_data$snp_id[order(snp_data$meanDP)])
ggplot(vcf, aes(x = SNPID, y = DP)) + geom_point() + xlab("SNP ID") + ylab("Depth")
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png)

SNP Mean Depth


```r
snp_data$SNPID <- factor(snp_data$snp_id, levels = snp_data$snp_id[order(snp_data$meanDP)])
ggplot(snp_data, aes(x = SNPID, y = meanDP)) + geom_point() + xlab("SNP ID") + ylab("Mean Depth")
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png)


Minor allele frequency histogram


```r
ggplot(snp_data, aes(x = maf)) + geom_histogram(bins = 20)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)

Genotype frequencies


```r
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(snp_data, aes(x = af)) + geom_point(aes(y = paa), color = cbPalette[1]) + geom_point(aes(y = pAa), color = cbPalette[2]) + geom_point(aes(y = pAA), color = cbPalette[3]) + xlab("Allele a frequency") + ylab("Genotype frequency")
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png)



Call rate x heterozygosity rate


```r
ggplot(snp_data, aes(x = pAa, y = call_rate)) + geom_point()
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png)

```r
ggplot(snp_data, aes(x = maf, y = call_rate)) + geom_point()
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-2.png)


Call rate x pcr product size


```r
ggplot(snp_data, aes(x = pcr1_product, y = call_rate)) + geom_point() + geom_smooth(method = "lm") + xlab("PCR 1 product size") + ylab("Call rate")
```

```
## `geom_smooth()` using formula 'y ~ x'
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png)


PCR product size x mean depth


```r
ggplot(snp_data, aes(x = pcr1_product, y = meanDP)) + geom_point() + ylab("Mean depth") + xlab("PCR 1 product size")
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png)


Call rate x mean depth


```r
ggplot(snp_data, aes(x = meanDP, y = call_rate)) + geom_point() + xlab("Mean depth") + ylab("Call rate")
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.png)


Call rate x sample type


```r
ggplot(vcf, aes(x = type, y = DP)) + geom_boxplot()
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png)



