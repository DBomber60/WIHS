setwd("~/Documents/personal/diss/WIHS")
library(tidyverse)

#drugmap = rbind(drugmap, c("", 998, "NRTI + NNRTI"))
#> drugmap = rbind(drugmap, c("", 999, "All"))

labs = read.csv("labsum.csv")
f22r = read.csv('f22r.csv') # f22r.csv: new recruit antiviral medication history
drug1 = read.csv('drug1.csv') # drug1.csv: antiviral medications
drugmap = read.csv('drugmap.csv', stringsAsFactors = FALSE) # drug code to classification

# subset to HIV + patients who have not taken an antiviral (EMTXDR=2)
# filter out HIV negative patients (maindr != 1)

f22r.f = f22r %>% filter(EMTXDR == 2, MAINDR != 1)
drug1.f = drug1 %>% filter(CASEID %in% f22r.f$CASEID)
# length(unique(drug1.f$CASEID)) # 104 (out of 132) end up being treated

# let's look at these treatment sequences/ how much data do we have here 
case.vis = drug1.f %>% group_by(CASEID, VISIT) %>% 
  summarise(n=n()) %>% mutate(start = min(VISIT)) %>% arrange(CASEID, VISIT)

# now, let's look at the kinds of drugs we observe in these sequences


which.drugs = inner_join(case.vis, drug1.f, by = c("CASEID", "VISIT")) %>%
  left_join(drugmap, by = c("DRUGD1"="Num")) %>% 
  dplyr::select(CASEID, VISIT, DRUGD1, Class)


l = which.drugs %>% group_by(CASEID, VISIT) %>% group_map(~ unique(.x$Class))

# lapply a function 

# l is a character vector with at least two elements
dosage = function(l) {
  
  if(length(l) == 1) {return(l)} else {
  
  if ( "NNRTI" %in% l & "NRTI" %in% l ) {
    return("NRTI + NNRTI")
  } else if ( "NNRTI" %in% l & "PI" %in% l ) {
    return("NNRTI + PI")
  } else if ( "NRTI" %in% l & "PI" %in% l ) {
    return("NRTI + PI")
  } else if ( "NRTI" %in% l & "II" %in% l ) {
    return("NRTI + II")
  } else {return("All")} }
}

case.vis$class = unlist( lapply(l, dosage) )
# map caseid to a number
case.vis = case.vis %>% filter(class != "NNRTI + PI", class != "All")
case.vis$id = sapply(case.vis$CASEID, function(x) which(unique(case.vis$CASEID)==x))


ggplot(case.vis, aes(VISIT, id)) + geom_raster(aes(fill=class)) + 
  theme_minimal() + ggtitle("Treatment Sequences") + ylab("Patient ID") + xlab("Visit Number")

#ggsave("treatseq.png", dpi = 320, width = 12, height = 7)

# TO DO (5.10)
# 1. make a stacked area chart (https://www.r-graph-gallery.com/stacked-area-graph.html)
# 2. add in the other untreated people to these TVC plots



#### area plot to show the proportion of each treatment
plot = case.vis %>% group_by(VISIT, class) %>% summarise(n = n()) %>% 
  mutate(prop = n/length(unique(case.vis$CASEID)))

ggplot(plot, aes(x = VISIT, y = prop, fill = class)) + geom_area() + theme_minimal() + 
  ggtitle("Antiretroviral Therapy Use Among WIHS HIV+ Patients") + 
  ylab("% Receiving Therapy") + xlab("Visit")

ggsave("therapy.png", dpi = 320, width = 12, height = 7)

# figure out how to get the untreated in here

# let's look at their lab values for this same time frame
labs.f = filter(labs, CASEID %in% unique(case.vis$CASEID))

# let's plot the time varying confounders for a single visit
labdat.plot = left_join(labs.f, case.vis, by = c("CASEID","VISIT")) %>% 
  filter(VISIT == 35)
labdat.plot$class = ifelse(is.na(labdat.plot$class), "Untreated", labdat.plot$class)

labdat.plot = filter(labdat.plot, class != "NRTI + II")

p1 = ggplot(labdat.plot, aes(x = class, y = log(VLOAD))) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_jitter(0.1)) + theme_minimal() + 
  ggtitle("HIV RNA (copies per ml)") + xlab("")

p2 = ggplot(labdat.plot, aes(x = class, y = log(CD4N))) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_jitter(0.1)) + theme_minimal() +
  ggtitle("Log(# of CD4 cells)") + xlab("")

p3 = ggplot(labdat.plot, aes(x = class, y = log(CD8N))) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_jitter(0.1)) + theme_minimal() + 
  ggtitle("Log(# of CD8 positive cells)") + xlab("")

p4 = ggplot(labdat.plot, aes(x = class, y = HSRAT)) + geom_boxplot() + 
  geom_jitter(shape=16, position=position_jitter(0.1)) + theme_minimal() + 
  ggtitle("Helper/ Suppressor Ratio") + xlab("")


g = arrangeGrob(p1, p2, p3, p4, nrow = 2)
ggsave("eda.pdf", g)

ggplot(labdat, aes(x=VISIT, y = CD4N, group = CASEID)) + geom_line(aes(color=CASEID))
ggplot(labdat, aes(x=VISIT, y = log(VLOAD), group = CASEID)) + geom_line(aes(color=CASEID))


# Box plot with dot plot
p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)
# Box plot with jittered points
# 0.2 : degree of jitter in x direction
p + geom_jitter(shape=16, position=position_jitter(0.2))


ds = which.drugs %>% group_by(DRUGD1) %>% summarise(n=n())

# 211 PI
# 253 NRTI
# 243 PI + Pharmacokinetic Enhancers
# 220 NNRTI
# 262 NRTI
# 227 NRTI
# 204 NRTI
# 234 NRTI
# 217 PI
# 191 NNRTI
# 240 NRTI
# 147 NRTI
# 254 NRTI
# 218 NRTI
# 216 PI
# 249 PI
# 280 NNRTI + NRTI
# 159 NRTI
# 239 NRTI
# 256 PI
# 264 II
# 265 EI
# 




# II: 286, 284, 264

# let's get labs data for these visits
case.labs = inner_join(seqdat2, labs, be = c("CASEID", "VISIT"))


ggplot(case.labs, aes(x=VISIT, y = CD4N, group = CASEID)) + 
  geom_line(aes(color=factor(CASEID)))

# delineate drug type by visit by patient

type.df = left_join(which.drugs, drugmap, by = c("DRUGD1" = "Num")) %>% select(CASEID, VISIT, Class)

# let's look at time of seroconversion as a time 0

hivdat = read.csv("hivstat.csv")

# people looked at before - were they positive or negative?
 
caseh = inner_join(seqdat2, hivdat, by = "CASEID") %>% group_by(CASEID) %>%
  summarise(m = max(STATUS), d = max(POSVIS), g = max(FPOSDATE))

# were these people actually not on treatment before? let's look at one example (1574)
# correct - per f22r, this person had never taken ART b/c her viral load was too low

# treatment types

# filter(f22r, CASEID == "1574")
t2 = filter(drug1, CASEID == "1574")

# next, let's not get rid of people for whom we do not have consecutive ART data.




