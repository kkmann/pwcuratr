### Prune pathway cluster

To increase the specificiy of the elicited pathway clusters, 
the candidate genes derived via shared reactome pathways with any of the
seed genes can be pruned.
The pruing is done on the graph defines by all candidate genes and the 
set of interacitons by Wu et al. filtered by a minimal confidence score.
Only candidate genes that are 
$k$ 
or fewer edges (interactions) away from one of
the seed genes are retained.
