# HPart
Greedy Randomized Adaptive Search Procedure for Network 2-partition

	HPART.C: A GRASP approach to the 2-partition problem

	GRASP == Greedy Randomized Adaptive Search Procedure
	Builds a low weight partition of 0-1 graph by greedily adding
	pairs of nodes from a candidate list of those nodes the maximize
	the current partition, then the weight of the partition is reduced
	by exchanging pairs of nodes when the exchange will increase the
	weight of the inner edges.

	Hal Elrod   --  Operations Research Group
			Department of Mechanical Engineering
			The University of Texas at Austin

	Spring Break, March 1989
	last modified 6/30/89
   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

   Input:	number of nodes, number of edges
		node, node, weight
		 .      .     .
		 .	.     .
		node, node, weight

   Use:		HPART <inputfile> <modea> <modeb> <c-list> <run-time>
   Where:	<modea> = "1", heap greedy unmatched partition
			= "2", greedy unmatched partition
		<modeb> = "1", first swap -- generic 2-exchange
			= "2", slight swap
			= "3", slightest swap
			= "4", compact slight swap
