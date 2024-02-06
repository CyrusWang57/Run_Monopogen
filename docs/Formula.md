# Issue about formula

## Methods in paper

The pictures below demonstrates the method of Monopogen described in paper.

<figure markdown>
![Image title](./monopogen_image/formula_monopogen_1.png){ width="80%" }
</figure>

<figure markdown>
![Image title](./monopogen_image/formula_monopogen_2.png){ width="80%" }
</figure>

In brief, for each custom-defined physical distance bin, Monopogen calculates a linkage proportion for germline SNV pairs within that distance, and then takes its reciprocal, meaning the closer this value is to 1, the more likely it is that the remaining allele of two SNVs is on the same haplotype within this distance range.  
Subsequently, this is used as a reference to calculate the phasing probability of the allele of a de novo SNV, which involves calculating the probability that it is on the same haplotype with adjacent germline SNVs for both phasing scenarios. An average value is taken across all cells, with the higher value indicating the final phasing at the cell population level, and the lower value serving as an indicator for measuring somatic mutations (the higher it is, the more likely it is somatic).

## Code implementation

The code below calculates the phasing probability of de novo SNV.  

```R
weigthedP <- function(dis=NULL, match=NULL, table=NULL){
  binsize <- c(0, 100, 200, 300, 400, 500, 1000, 2500, 5000, 7500, 10000, 20000, 50000, 100000, 500000, 1000000000000000000)
  rd_p <-c()
  rd_w <-c()
  for(i in seq(2,length(binsize),1)){
    index <-which(dis<binsize[i] & dis>=binsize[i-1])
    if(length(index>0)){
      pval <- sum(match[index])/length(index)
      rd_p <- c(rd_p, pval)
      rd_w <- c(rd_w, (1-table$Prob[i-1])^2)
    }
  }
  p <- 1-sum(rd_p*rd_w)/(sum(rd_w))
  return(p)
}
```


\begin{align}
M_{b_i} = \frac {\sum_{c_i \in C} (V_{match} | b_i, c_i)} {N_C} \\
P_{b_i} = (1 - P_{two_loci} | b_i)^2 \\
P_{final} = 1 - \frac {\sum_{b_i \in B} (M_{b_i} * P_{b_i})} {\sum_{b_i \in B} P_{b_i}}
\end{align}

The homomorphism $f$ is injective if and only if its kernel is only the
singleton set $e_G$, because otherwise $\exists a,b\in G$ with $a\neq b$ such
that $f(a)=f(b)$.