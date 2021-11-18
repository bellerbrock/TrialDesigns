Trial Design Options
================

#### Table of Contents

[Randomized Complete Block
Design](#randomized-complete-block-design)<br> [Resolvable Row Column
Design](#resolvable-row-column-design)<br> [Doubly Resolvable Row Column
Design](#doubly-resolvabl-row-column-design)<br> [Field
Dimensions](#field-dimensions)<br> [Field
Variation](#field-variation)<br>

## Randomized Complete Block Design

Typical 96 plot (24 accession, 4 rep) randomized complete block design
(CB) using agricolae package. No restriction on randomization within the
one-way blocks. Simple, but can result in significant clumping.

``` r
library(agricolae)
library(desplot)

checks <- c("Covington", "Beauregard")
treatments <- c(checks, 1:22)
nRep <- 4
nRow <- 8
nTier <- 12
plotWidth <- 7  # 7 foot wide plots (42 inch field rows x 2 field rows/plot)
plotLength <- 25
xlen <- nRow * plotWidth
ylen <- nTier * plotLength

CB <- design.rcbd(treatments, r=nRep, serie=2)$book

CB <- transform(CB, is_a_control = ifelse(CB$treatments %in% checks, TRUE, FALSE))
CB$row = rep(1:nRow,times=nTier)
CB$tier = rep(1:nTier,each=nRow)

# ggplot display, excludes treatments from key by default
CB_gg_layout <- desplot( form=block~row+tier, data=CB, text=treatments, out1=block,
  col=is_a_control, cex=1.25, main="Randomized Complete Block Design",
  gg=TRUE, xlab="rows", ylab="tiers", shorten="sub")
CB_gg_layout
```

![](README_files/figure-gfm/agricolae%20CB-1.png)<!-- -->

## Resolvable Row Column Design

96 plot (24 accession, 4 rep) resolvable row-column design (RC) using
blocksdesign package. In addition to the one-way blocks it adds a
constraint that no treatment may appear more than once in the long rows
(latinization). Each long row becomes an incomplete block.

``` r
library(blocksdesign)

RCblocks <- data.frame(
  block = gl(nRep,length(treatments)),
  row = gl(nRow,1),
  tier = gl(nTier,nRow)
)

RC <- design(treatments, RCblocks)$Design
RC <- transform(RC, is_a_control = ifelse(RC$treatments %in% checks, TRUE, FALSE))
RClayout <- desplot(form=block~row+tier, data=RC, text=treatments, out1=block,
    out2=row, col=is_a_control, cex=1.25, main="Resolvable Row-Column Design",
    gg=TRUE, xlab="rows", ylab="tiers", shorten="sub")
RClayout
```

![](README_files/figure-gfm/blocksdesign%20RC-1.png)<!-- -->

## Doubly Resolvable Row Column Design

96 plot (24 accession, 4 rep) doubly resolvable row-column design (DR)
using blocksdesign package. Like the resolvable row-column design, but
adds a constraint that no treatment may appear more than once in each
group of t=nRows/nReps long rows (t-latinization). Each group of long
rows forms a complete block, so the design has 2-way blocking and is
doubly resolvable.

``` r
create_dr_layout <- function(treatments, checks, reps, rows, tiers) {
  DRblocks <- data.frame(
    tier_block = gl(reps,length(treatments)),
    row_block = gl(reps,1),
    dummy_tier = gl(tiers*2,rows/2)
  )

  DR <- design(treatments, DRblocks)$Design
  
  DR <- transform(DR, is_a_control = ifelse(DR$treatments %in% checks, TRUE, FALSE))
  DR$row = rep(rep(1:rows,each=tiers/reps),times=reps)
  DR <- DR[order(DR$row),]
  DR$tier = rep(1:tiers,times=rows)
  
  return(DR)
}

DR <- create_dr_layout(treatments, checks, reps=nRep, rows=nRow, tiers=nTier)

DRlayout <- desplot(form=tier_block~row+tier, data=DR, text=treatments, out1=tier_block,
    out2=row_block, col=is_a_control, cex=1.25, main="Doubly Resolvable Row-Column Design",
    gg=TRUE, xlab="rows", ylab="tiers", shorten="sub")
DRlayout
```

![](README_files/figure-gfm/blocksdesign%20DR-1.png)<!-- -->

## Field Dimensions

Comparison of field dimensions for different row and tier configurations
using desplot aspect option. Plots are longer than they are wide, so a
greater number of rows than tiers should make squarer designs, and be
more suitable for 2-way blocking.

``` r
plotWidth <- 7  # 42 inch spacing between field rows, 2 field rows per plot = 7 foot wide plots
plotLength <- 25
xlen <- nRow * plotWidth
ylen <- nTier * plotLength

# lattice display instead of gg, incorporates field dimensions via aspect=ylen/xlen
lattice_layout <-desplot(form=tier_block~row+tier, data=DR, text=treatments, out1=tier_block,
    out2=row_block, col=is_a_control, main=paste(nRow," rows X ", nTier," tiers"),
    aspect=ylen/xlen, xlab="rows", ylab="tiers", shorten="sub", show.key=FALSE, cex=.65)

# Switch number of rows and tiers and compare
altDR <- create_dr_layout(treatments, checks, reps=nRep, rows=nTier, tiers=nRow)
altxlen <- nTier * plotWidth
altylen <- nRow * plotLength
lenratio = nRow/nTier

alt_lattice_layout <- desplot(form=tier_block~row+tier, data=altDR, text=treatments,
    out1=tier_block, out2=row_block, col=is_a_control, main=paste(nTier," rows X ", nRow," tiers"),
    aspect=altylen/altxlen, xlab="rows", ylab="tiers", shorten="sub", show.key=FALSE, cex=.65)

library(cowplot)
plot_grid(lattice_layout, alt_lattice_layout, scale = c(1, lenratio))
```

<img src="README_files/figure-gfm/desplot aspect-1.png" style="display: block; margin: auto;" />

## Field Variation

``` r
# library(brapir)
# ??brapir
# spbase <- brapi_db()$sweetpotatobase
# T17PRE10024HCR = brapi_get_studies_studyDbId(con = spbase, studyDbId = "1048")
# T17PRE10024HCR <- brapi_get_studies_studyDbId_observationunits(con = spbase, studyDbId = "1048")

T18NCGT0014HCR <- read.csv("18NCGT0014HCR.csv", skip =3, fill=TRUE, sep=",",
        header = TRUE, stringsAsFactors = T, na.strings="NA");
# colnames(T18NCGT0014HCR)
names(T18NCGT0014HCR)[names(T18NCGT0014HCR) ==
        "Total.storage.root.weight.per.NET.plot.in.kg.CO_331.0000237"] <- "TRW"
desplot(T18NCGT0014HCR, TRW ~ colNumber*rowNumber, main="18NCGT0014HCR Total Root Yield",
        col=germplasmName, num=germplasmName, cex=1, out1=blockNumber, aspect=4/1)
```

<img src="README_files/figure-gfm/desplot heatmap-1.png" style="display: block; margin: auto;" />

``` r
# names(T18NCGT0014HCR)[names(T18NCGT0014HCR) ==
# "Length.to.diameter.ratio.computation.CO_331.0000779"] <- "LDR"
# desplot(T18NCGT0014HCR, LDR ~ colNumber*rowNumber, main="18NCGT0014HCR Length to Diameter Ratio",
#         col=germplasmName, num=germplasmName, cex=1, out1=blockNumber, aspect=4/1)
```
