##
# Port Excel power calculations to R from Excel
##
library(pwr)
library(dplyr)
library(ggplot2)
library(tidyr)
#
A = 1 ## Additive Effect
D = 1 ## Dominance Effect
p = 0.025 # Minor Allele frequency
h2 = 0.249 # heritability (excluding this locus)
f = 0.001 # proportion of phenotypic variance explained by this locus
g = 100 # Number of Daughters per sire
s = 5e-8 # Desired Significance (Power level)
x = 1 - s # for critical value calculations
#
# genotype frequencies(gf) (rr, rm, mm; reference-reference, reference-mutant, mutant-mutant)
gf.rr = (1-p)^2
gf.rm = 2*p*(1-p)
gf.mm = p^2
#
# genotypic values(gv)
gv.rr = A
gv.mm = -A
gv.rm = ((gv.rr + gv.mm)/2) + D
# Trait values
t.mean = gf.rr*gv.rr + gf.rm*gv.rm + gf.mm*gv.mm # trait mean
t.vl = gf.rr*(gv.rr-t.mean)^2 + gf.rm*(gv.rm-t.mean)^2 + gf.mm*(gv.mm-t.mean)^2 # trait variance locus
t.vo = (1-f)*t.vl/f # trait variance other
t.vg = h2*t.vo # trait variance genetic other
t.ve = (1-h2)*t.vo # trait variance environment
t.vy = t.vo + t.vl # phenotypic variance
# allelic effects(ae) (r - reference, m - mutant)
ae.r = (1-p)*gv.rr + p*gv.rm - t.mean
ae.m = (1-p)*gv.rm + p*gv.mm - t.mean
# breeding values(bv)
bv.rr = 2*ae.r + t.mean
bv.rm = ae.r + ae.m + t.mean
bv.mm = 2*ae.m + t.mean
# Power function, parameters = # Samples, # Model Parameters, Phenotypic variance, Significance threshold
get_power = function(n, k, t.vy, vl.m, s)
{
  return(pwr.f2.test(k-1, n-k, vl.m/(t.vy), sig.level = s)$power)
}
#
dataset = tibble(N = seq(500, 100000, 500))
# Cows (cow)
# additive model (a)
cow.a.vl.m = gf.rr*(bv.rr-t.mean)^2 + gf.rm*(bv.rm-t.mean)^2 + gf.mm*(bv.mm-t.mean)^2 # locus variance in the model (vl.m)
cow.a.vl.n = t.vl - cow.a.vl.m # locus variance not in the model (vl.n)
cow.a.vr = t.vg + t.ve + cow.a.vl.n # residual variance (vr)
cow.a.k = 2 # number of parameters in the model (k)
dataset = dataset %>% mutate(cow.a.power = get_power(N, cow.a.k, t.vy, cow.a.vl.m, s)) # Power via pwr function instead
# recessive model (r)
cow.r.vl.m = (gf.rr+gf.rm)*(((gf.rr*gv.rr + gf.rm*gv.rm)/(gf.rr+gf.rm)) - t.mean)^2 + gf.mm*(gv.mm - t.mean)^2
cow.r.vl.n = t.vl - cow.r.vl.m
cow.r.vr = t.vg + t.ve + cow.r.vl.n
cow.r.k = 2
dataset = dataset %>% mutate(cow.r.power = get_power(N, cow.r.k, t.vy, cow.r.vl.m, s))
# genotype class effect model (g)
cow.g.vl.m = gf.rr*(gv.rr-t.mean)^2 + gf.rm*(gv.rm-t.mean)^2 + gf.mm*(gv.mm-t.mean)^2
cow.g.vl.n = t.vl - cow.g.vl.m
cow.g.vr = t.vg + t.ve + cow.g.vl.n
cow.g.k = 3
dataset = dataset %>% mutate(cow.g.power = get_power(N, cow.g.k, t.vy, cow.g.vl.m, s))
#
# Bulls (bull)
# additive model (a)
bull.a.vl.m = gf.rr*(ae.r^2) + gf.rm*((0.5*(ae.r+ae.m))^2) + gf.mm*(ae.m^2) 
bull.a.vl.n = 0 # the additive model explains everything for the locus
bull.a.vr = (0.75*t.vg+t.ve)/g
bull.vy = bull.a.vr + bull.a.vl.m # bull phenotypic variance
bull.a.k = 2
dataset = dataset %>% mutate(bull.a.power = get_power(N, bull.a.k, bull.vy, bull.a.vl.m, s))
# recessive model (r)
bull.r.vl.m = (gf.rm+gf.rr)*((gf.rr*ae.r + gf.rm*0.5*(ae.r+ae.m))/(gf.rm+gf.rr))^2 + gf.mm*(ae.m^2)
bull.r.vl.n = bull.a.vl.m - bull.r.vl.m
bull.r.vr = (0.75*t.vg+t.ve)/g + bull.r.vl.n
bull.r.k = 2
dataset = dataset %>% mutate(bull.r.power = get_power(N, bull.r.k, bull.vy, bull.r.vl.m, s))
# genotype class effect model (g)
bull.g.vl.m  = gf.rr*(ae.r^2) + gf.rm*((0.5*(ae.r+ae.m))^2) + gf.mm*(ae.m^2)
bull.g.k = 3
dataset = dataset %>% mutate(bull.g.power = get_power(N, bull.g.k, bull.vy, bull.g.vl.m, s))
# additive model - no homozygous alternates (q)
bull.q.mean = (gf.rr*(t.mean + bv.rr/2) + gf.rm*(t.mean + bv.rm/2))/(1-gf.mm)
bull.q.vl.m = (gf.rr*(t.mean - bull.q.mean + bv.rr/2)^2 + gf.rm*(t.mean - bull.q.mean + bv.rm/2)^2)/(1-gf.mm)
bull.q.vl.n = 0
bull.q.vr = (0.75*t.vg+t.ve)/g + bull.q.vl.n
bull.q.k = 2
dataset = dataset %>% mutate(bull.q.power = get_power(N, bull.q.k, bull.vy, bull.q.vl.m, s))
#
dataset.power = select(dataset, N, ends_with("power"))
colnames(dataset.power) = c("N", "Cow - Additive", "Cow - Recessive", "Cow - Genotype Class", "Bull - Additive", "Bull - Recessive", "Bull - Genotype Class", "Bull - Additive - No affecteds")
dataset.power.g = gather(dataset.power, model, Power, -N)
dataset.power.g = mutate(dataset.power.g, model = factor(model, levels = c("Cow - Additive", "Cow - Recessive", "Cow - Genotype Class", "Bull - Additive", "Bull - Recessive", "Bull - Genotype Class", "Bull - Additive - No affecteds")))
#
ggplot(dataset.power.g, aes(N, Power, col = model)) + geom_line()

