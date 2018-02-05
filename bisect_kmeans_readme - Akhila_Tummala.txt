BISECT K-MEANS:
	This program has the code for the bisect K-Means algorithm. It randomly generates the specified number of points in the specified dimensions.
	Then, the points are assigned to the clusters based on the bisecting k-means algorithm used for calculating the initial cluster centroids and
	then the k-means algorithm is applied on these initial clusters.
	Then, 10 random points are generated and these 10 points, one by one, are searched for their neighbors in these clusters.
	The output of this program contains:
			1. Minimum distance between each random point and its nearest point
			2. Number of points traversed in order to find the minimum distance

Compile Instructions:
	The input C file takes three integer arguments as input:
		1. No. of dimensions
		2. No. of points
		3. No. of clusters
	The C file can be compiled using the following command:
		$gcc "bisect_kmeans - Shashank_Cheruku_Akhila_Tummala.c" -o bisect
	The output can be looked by executing the following command:
		$./bisect

Sample Input:
	The sample input is given while running the code.
	$./bisect
		Enter the number of dimensions: 8
		Enter the number of points: 1000
		Enter the number of clusters: 12
	
Sample Output:
	KD Tree search results:
        For random point 1:
                Total Points Visited: 4
                Minimum Distance: 25.573424
        For random point 2:
                Total Points Visited: 30
                Minimum Distance: 18.734994
        For random point 3:
                Total Points Visited: 20
                Minimum Distance: 14.071247
        For random point 4:
                Total Points Visited: 22
                Minimum Distance: 24.207437
        For random point 5:
                Total Points Visited: 19
                Minimum Distance: 21.679483
        For random point 6:
                Total Points Visited: 14
                Minimum Distance: 16.401219
        For random point 7:
                Total Points Visited: 6
                Minimum Distance: 14.628739
        For random point 8:
                Total Points Visited: 54
                Minimum Distance: 23.958297
        For random point 9:
                Total Points Visited: 8
                Minimum Distance: 30.626786
        For random point 10:
                Total Points Visited: 16
                Minimum Distance: 19.235384