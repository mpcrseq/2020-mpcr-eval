
# mpcrseq exploratory

Import data and join amplicon size data with genotyping data.


```r
library(tidyverse)
library(lme4)

vcf_129 <- read_csv("vcf_129snp_long.csv")
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
vcf_av  <- read_tsv("vcf_av_snps_long.txt")
```

```
## Parsed with column specification:
## cols(
##   CHROM = col_double(),
##   POS = col_double(),
##   REF = col_character(),
##   ALT = col_character(),
##   QUAL = col_double(),
##   FORMAT = col_character(),
##   sample = col_character(),
##   GT = col_character(),
##   PL = col_character(),
##   DP = col_double(),
##   AD = col_character(),
##   GQ = col_double()
## )
```

```r
vcf_129 <- vcf_129 %>% select(CHROM, POS, REF, ALT, sample = SAMPLE_NAME_IN_RUN, GT, PL, DP, AD, GQ)
vcf_av  <- vcf_av %>% select(CHROM, POS, REF, ALT, sample, GT, PL, DP, AD, GQ)
vcf     <- rbind(vcf_129, vcf_av)

# Separate allele depth & genotype likelihoods
vcf <- vcf %>% separate(GT, c("GT1", "GT2"), "/", convert = TRUE, remove = FALSE) %>%
  separate(AD, c("AD1", "AD2"), ",", convert = TRUE, remove = FALSE) %>%
  separate(PL, c("pl00", "pl01", "pl11"), convert = TRUE, remove = FALSE) %>%
  rowwise() %>% mutate(gt = sum(GT1, GT2))
```

```
## Warning: Expected 2 pieces. Missing pieces filled with `NA` in 3156 rows [32, 48, 72, 81, 124, 161, 177, 201, 210, 253, 290, 306, 330, 339, 382, 419, 435, 459, 468, 511, ...].
```

```
## Warning: Expected 3 pieces. Missing pieces filled with `NA` in 3156 rows [32, 48, 72, 81, 124, 161, 177, 201, 210, 253, 290, 306, 330, 339, 382, 419, 435, 459, 468, 511, ...].
```

```r
# Calc. likelihood ratio
likelihood_ratio <- function(x) {
  y <- 10^(-x/10)
  max(y) / (sum(y) - max(y))
}

vcf <- vcf %>% rowwise() %>% mutate(lr = likelihood_ratio(c(pl00,pl01,pl11)))




load("run0802.RData")
load("run1012.RData")

r1 <- run0802$run_x_sample_snp_data
r2 <- run1012$run_x_sample_snp_data

r12 <- rbind(r1, r2)

r12 <- r12 %>% select(CHROM = chrom, POS = stop, REF = ref, ALT = alt, sample, gt, pl00, pl01, pl11, DP = dp, gt, lr)

vcf <- vcf %>% select(CHROM, POS, REF, ALT, sample, gt, pl00, pl01, pl11, DP, gt, lr)


vcf <- rbind(r12, vcf)

s1s2        <- read_tsv("set12_product_sizes.tsv")
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
vcf        <- left_join(s1s2, vcf)
```

```
## Joining, by = "snp_id"
```

```r
vcf        <- vcf %>% filter(!is.na(CHROM))

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
vcf      <- left_join(vcf, metadata)
```

```
## Joining, by = "sample"
```

```r
pchar       <- read_tsv('pchar.tsv')
```

```
## Parsed with column specification:
## cols(
##   chrom = col_double(),
##   start = col_double(),
##   stop = col_double(),
##   amplicon_snp_position = col_double(),
##   l_start = col_double(),
##   r_start = col_double(),
##   l_length = col_double(),
##   r_length = col_double(),
##   amplicon_size = col_double(),
##   f_seq = col_character(),
##   r_seq = col_character(),
##   amplicon = col_character(),
##   template = col_character(),
##   f_gc_clamp = col_logical(),
##   r_gc_clamp = col_logical(),
##   f_gc_3prime = col_logical(),
##   r_gc_3prime = col_logical()
## )
```

```r
pchar$start <- as.integer(pchar$start)
pchar$name  <- paste(pchar$chrom, pchar$start, sep = ':')
pchar       <- pchar %>% select(snp_id = name, f_gc_clamp, r_gc_clamp, f_gc_3prime, r_gc_3prime, f_seq, r_seq) %>% distinct()
pchar       <- pchar %>% mutate(gc_clamp = f_gc_clamp & r_gc_clamp, gc_3prime = f_gc_3prime, r_gc_3prime)
vcf         <- left_join(vcf, pchar)
```

```
## Joining, by = "snp_id"
```

Generate some variables, and group data by snp_id.


```r
vcf$called <- vcf$lr > 10
# Group by snp_id, calc snp-level variables.

snp_data <- vcf %>% 
  group_by(snp_id, CHROM, REF, ALT, template_size, pcr1_product, pcr2_product, gc_clamp, gc_3prime, f_gc_clamp, r_gc_clamp, f_gc_3prime, r_gc_3prime) %>% 
  summarize(n_called = sum(called), n = n(), call_rate = n_called/n, totalDP = sum(DP), meanDP = mean(DP), medianDP = median(DP),
  paa = sum(gt == 0) / n(), pAa = sum(gt == 1) / n(), pAA = sum(gt == 2) / n(), a1 = (2*sum(!is.na(gt))) - sum(gt, na.rm = TRUE), a2 = sum(gt, na.rm = TRUE))

multiallelic <- snp_data %>% group_by(snp_id) %>% filter(n() > 1)

snp_data <- snp_data %>% filter(!snp_id %in% multiallelic$snp_id)

# Minor allele frequency

getmaf <- function(a1, a2) {
  if(a1 == 0 || is.na(a1) || a2 == 0 || is.na(a2)) {return(0)}

  af1 <- a1 / (a1 + a2)
  if(af1 > 0.5){
    return(1 - af1)
  } else {
    return(af1)
  }
}

snp_data <- snp_data %>% rowwise() %>% mutate(maf = getmaf(a1, a2))

# Allele frequency
getaf <- function(a1, a2) {

  if(a2 == 0 || is.na(a2)) {return(0)}
  if(a1 == 0 || is.na(a1)) {return(1)}

  af1 <- a1 / (a1 + a2)
}

snp_data <- snp_data %>% rowwise() %>% mutate(af = getaf(a1, a2))

sample_data <- vcf %>% group_by(sample, type, dna, Year, region, run_id) %>% summarise(n_called = sum(called), n = n(), call_rate = n_called/n, totalDP = sum(DP), meanDP = mean(DP), medianDP = median(DP))

sample_data <- sample_data %>% group_by(run_id) %>% mutate(totalRunDP = sum(totalDP)) %>% mutate(scaledDP = totalDP / totalRunDP)
```
## Figure 1. SNP Depth



```r
vcf$SNPID <- factor(vcf$snp_id, levels = snp_data$snp_id[order(snp_data$meanDP)])
ggplot(vcf, aes(x = SNPID, y = DP)) + geom_point() + xlab("SNP ID") + ylab("Depth") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

```
## Warning: Removed 1289 rows containing missing values (geom_point).
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png)

## Figure 2. SNP Mean Depth


```r
snp_data$SNPID <- factor(snp_data$snp_id, levels = snp_data$snp_id[order(snp_data$meanDP)])
ggplot(snp_data, aes(x = SNPID, y = meanDP)) + geom_point() + xlab("SNP ID") + ylab("Mean Depth") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

```
## Warning: Removed 4 rows containing missing values (geom_point).
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png)

## Figure 3. Minor allele frequency histogram


```r
ggplot(snp_data, aes(x = maf)) + geom_histogram(bins = 20)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)

## Figure 4. Genotype frequencies


```r
p <- seq(0,1,length.out = nrow(snp_data))
q <- 1-p
pq = (1-(p^2)-(q^2))/2

aa = p^2
AA = q^2
Aa = 2*pq

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(snp_data, aes(x = af)) + 
  geom_point(aes(y = paa), color = cbPalette[1]) + 
	geom_line(aes(x = p, y = aa), color = cbPalette[1]) +
	geom_point(aes(y = pAa), color = cbPalette[2]) + 
	geom_line(aes(x = p, y = Aa), color = cbPalette[2]) +
	geom_point(aes(y = pAA), color = cbPalette[3]) + 
	geom_line(aes(x = p, y = AA), color = cbPalette[3]) +
	xlab("Allele a frequency") + ylab("Genotype frequency")
```

```
## Warning: Removed 178 rows containing missing values (geom_point).

## Warning: Removed 178 rows containing missing values (geom_point).

## Warning: Removed 178 rows containing missing values (geom_point).
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png)

## Figure 5. Call rate x heterozygosity rate


```r
ggplot(snp_data, aes(x = pAa, y = call_rate)) + geom_point()
```

```
## Warning: Removed 178 rows containing missing values (geom_point).
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png)

```r
ggplot(snp_data, aes(x = maf, y = call_rate)) + geom_point()
```

```
## Warning: Removed 44 rows containing missing values (geom_point).
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-2.png)

## Figure 6. Call rate x pcr product size


```r
ggplot(snp_data, aes(x = pcr1_product, y = call_rate)) + geom_point() + geom_smooth(method = "lm") + xlab("PCR 1 product size") + ylab("Call rate")
```

```
## `geom_smooth()` using formula 'y ~ x'
```

```
## Warning: Removed 44 rows containing non-finite values (stat_smooth).
```

```
## Warning: Removed 44 rows containing missing values (geom_point).
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png)

```r
glm1 <- glmer(called ~ I(scale(pcr1_product)) + (1|sample), family = binomial, data = vcf)
summary(glm1)
```

```
## Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
##  Family: binomial  ( logit )
## Formula: called ~ I(scale(pcr1_product)) + (1 | sample)
##    Data: vcf
## 
##      AIC      BIC   logLik deviance df.resid 
## 115397.5 115426.9 -57695.8 115391.5   134222 
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -13.6568  -0.4487  -0.0798   0.4855  12.4509 
## 
## Random effects:
##  Groups Name        Variance Std.Dev.
##  sample (Intercept) 5.062    2.25    
## Number of obs: 134225, groups:  sample, 1049
## 
## Fixed effects:
##                        Estimate Std. Error z value Pr(>|z|)    
## (Intercept)            -0.26864    0.07005  -3.835 0.000126 ***
## I(scale(pcr1_product)) -0.69152    0.00817 -84.641  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr)
## I(scl(p1_)) 0.005
```


## Figure 7. PCR product size x Mean depth


```r
ggplot(snp_data, aes(x = pcr1_product, y = meanDP)) + geom_point() + ylab("Mean depth") + xlab("PCR 1 product size")
```

```
## Warning: Removed 4 rows containing missing values (geom_point).
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png)


## Figure 8. Mean depth x Call rate


```r
ggplot(snp_data, aes(x = meanDP, y = call_rate)) + geom_point() + xlab("Mean depth") + ylab("Call rate")
```

```
## Warning: Removed 44 rows containing missing values (geom_point).
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.png)


## Figure 9. Sample type x Call rate


```r
ggplot(sample_data, aes(x = type, y = call_rate)) + geom_boxplot()
```

```
## Warning: Removed 1049 rows containing non-finite values (stat_boxplot).
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png)


## Figure 10 a,b,c. GC Clamp x Call rate

GC clamp = TRUE if 2 or fewer GC bases in 3 prime end.


>    Returns True if fewer than 2 G or C nucleotides are within the last 5 bases.
>    Oligonucleotides with more than three G or C nucleotides within the last  five bases were  filtered.  This should help minimize mispriming at GC-rich binding sites
>    http://www.bio.net/mm/methods/2000-August/084438.html



```r
ggplot(snp_data, aes(x = gc_clamp, y = call_rate)) + geom_boxplot() + ggtitle("GC Clamp")
```

```
## Warning: Removed 44 rows containing non-finite values (stat_boxplot).
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-1.png)

```r
ggplot(snp_data, aes(x = f_gc_clamp, y = call_rate)) + geom_boxplot() + ggtitle("Forward GC Clamp")
```

```
## Warning: Removed 44 rows containing non-finite values (stat_boxplot).
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-2.png)

```r
ggplot(snp_data, aes(x = r_gc_clamp, y = call_rate)) + geom_boxplot() + ggtitle("Reverse GC Clamp")
```

```
## Warning: Removed 44 rows containing non-finite values (stat_boxplot).
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-3.png)


## Figure 11 a,b,c. GC 3prime x Call rate

TRUE if GC is last base at 3prime end.

>    Rationale:
>    For primers, a single G or C nucleotide at the 3'-end helps to stabilize
>    binding near the site of extension, which can reduce the possibility of
>    "breathing" and improves priming efficiency.  terefore, primers ending
>    in an A or a T base could optionally be  altered. In this study, the 
>    filter was applied. [Francis et al. 2017 - ThermoAlign: a genome-aware
>    primer design tool for tiled amplicon resequencing]
>    http://www.bio.net/mm/methods/2000-August/084438.html



```r
ggplot(snp_data, aes(x = gc_3prime, y = call_rate)) + geom_boxplot() + ggtitle("GC 3prime")
```

```
## Warning: Removed 44 rows containing non-finite values (stat_boxplot).
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13-1.png)

```r
ggplot(snp_data, aes(x = f_gc_3prime, y = call_rate)) + geom_boxplot() + ggtitle("Forward GC 3prime")
```

```
## Warning: Removed 44 rows containing non-finite values (stat_boxplot).
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13-2.png)

```r
ggplot(snp_data, aes(x = r_gc_3prime, y = call_rate)) + geom_boxplot() + ggtitle("Reverse GC 3prime")
```

```
## Warning: Removed 44 rows containing non-finite values (stat_boxplot).
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13-3.png)


```r
glm8 <- glmer(called ~ f_gc_clamp + r_gc_clamp + f_gc_3prime + r_gc_3prime + (1|sample) + (1|snp_id), data = vcf, family = binomial )
summary(glm8)
```

```
## Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
##  Family: binomial  ( logit )
## Formula: called ~ f_gc_clamp + r_gc_clamp + f_gc_3prime + r_gc_3prime +      (1 | sample) + (1 | snp_id)
##    Data: vcf
## 
##      AIC      BIC   logLik deviance df.resid 
##  21463.1  21522.9 -10724.5  21449.1    37793 
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -223.756   -0.220   -0.013    0.221   41.071 
## 
## Random effects:
##  Groups Name        Variance Std.Dev.
##  sample (Intercept) 12.393   3.520   
##  snp_id (Intercept)  5.968   2.443   
## Number of obs: 37800, groups:  sample, 300; snp_id, 126
## 
## Fixed effects:
##                 Estimate Std. Error z value Pr(>|z|)  
## (Intercept)       0.7306     0.5498   1.329   0.1839  
## f_gc_clampTRUE   -1.1074     0.4600  -2.408   0.0161 *
## r_gc_clampTRUE   -0.8040     0.4911  -1.637   0.1016  
## f_gc_3primeTRUE  -0.6772     0.4704  -1.439   0.1500  
## r_gc_3primeTRUE   0.4481     0.4935   0.908   0.3639  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) f__TRU r__TRU f__3TR
## f_gc_clTRUE -0.448                     
## r_gc_clTRUE -0.485 -0.118              
## f_gc_3pTRUE -0.478  0.277  0.035       
## r_gc_3pTRUE -0.518 -0.030  0.419  0.069
```



