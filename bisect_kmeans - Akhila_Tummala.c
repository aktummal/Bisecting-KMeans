#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<limits.h>

void bisect_kmeans(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double *cluster_radius, double **cluster_centroid, int *cluster_assign);
void kmeans(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double *cluster_radius, double **cluster_centroid, int *cluster_assign);
void search_kmeans(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double *cluster_radius, double **cluster_centroid, double *query);


void main(){
	int dim=0, ndata=0, k=0, i=0, j=0, x=0, y=0, index=-1;
	double sum = 0.0, max = -1.0, min = -1.0;
	printf("Enter the number of dimensions: ");
	scanf("%d",&dim); // Dimensions input
	printf("Enter the number of points: ");
	scanf("%d",&ndata); //No. of points input
	printf("Enter the number of clusters: ");
	scanf("%d",&k); //No. of clusters input
	double *data; //Data array to store the points
	data = (double*) malloc((dim*ndata) * sizeof(double));
	int *cluster_size; //To store the sizes of the clusters
	cluster_size = (int*) malloc(k * sizeof(int));
	int *cluster_start; //To store the starting indexes of the clusters
	cluster_start = (int*) malloc(k * sizeof(int));
	double *cluster_radius; //To store the radius of the clusters
	cluster_radius = (double*) malloc(k * sizeof(double));
	int *cluster_assign; //To store the cluster number of each data point
	cluster_assign = (int*) malloc(ndata * sizeof(int));
	double **cluster_centroid; //To store the centroid of each cluster
	cluster_centroid = (double**) malloc(k * sizeof(double*));
	for(i=0;i<k;i++)
		cluster_centroid[i] = (double*) malloc(dim * sizeof(double));
	for(i=0;i<k;i++){ //Initializing the values to 0 initially
		cluster_size[i]=0;
		cluster_start[i]=0;
		cluster_radius[i]=0.0;
	}
	for(i=0;i<ndata;i++) //Initializing the values to -1
		cluster_assign[i]=-1;
	for(i=0;i<k;i++) //Initializing the values to 0
		for(j=0;j<dim;j++)
			cluster_centroid[i][j]=0.0;
	srand(time(NULL)); 
	for(i=0;i<dim*ndata;i++) // generate data set
		data[i] = rand() % 101;
	cluster_size[0] = ndata; //Initial size of the only cluster i.e. cluster 0
	bisect_kmeans(dim, ndata, data, k, cluster_size, cluster_start, cluster_radius, cluster_centroid, cluster_assign);
}

void bisect_kmeans(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double *cluster_radius, double **cluster_centroid, int *cluster_assign){
	int rand_index=-1, far_index=-1, total_clusters=1, i=0, j=0, parent=0, brk_cond=0, x=0;
	
	//<-------- This loop executes until the required no. of clusters are formed ------->
	while(total_clusters!=k){
		srand(time(NULL));
		rand_index = rand() % cluster_size[parent]; // Choosing the random point
		for(i=0;i<dim;i++)
			cluster_centroid[parent][i] = data[cluster_start[parent]+rand_index+i]; //Assigning the random point as cluster centroid
		double max_dist = LONG_MIN;
		
		//<----- Calculating the farthest point from the random point in that cluster ---->
		for(i=0;i<cluster_size[parent];i++){
			double temp_dist = 0.0;
			for(j=0;j<dim;j++)
				temp_dist = temp_dist + pow(data[cluster_start[parent]+i*dim+j]-cluster_centroid[parent][j],2);
			temp_dist = sqrt(temp_dist);
			if(temp_dist>max_dist){
				max_dist = temp_dist;
				far_index = cluster_start[parent] + i;
			}
		}
		
		for(i=0;i<dim;i++)
			cluster_centroid[total_clusters][i] = data[far_index+i]; //Assigning the farthest point as the other cluster centroid
		int temp_assign[cluster_size[parent]]; // Store the current current assigns before calculating the new ones
		for(i=0;i<cluster_size[parent];i++)
				temp_assign[i] = -1;
			
		//<------- Infinite loop breaks when the current and previous cluster centroids are equal ------>	
		while(1){
			
			//<----- Calculating the distances of each point from both cluster centroids and assigning the point to nearest cluster ----->
			for(i=0;i<cluster_size[parent];i++){
				double dist1=0.0, dist2=0.0;
				for(j=0;j<dim;j++){
					dist1 = dist1 + pow(data[cluster_start[parent]+i*dim+j]-cluster_centroid[parent][j],2);
					dist2 = dist2 + pow(data[cluster_start[parent]+i*dim+j]-cluster_centroid[total_clusters][j],2);
				}
				if(dist1<dist2)
					cluster_assign[cluster_start[parent]+i] = parent;
				else
					cluster_assign[cluster_start[parent]+i] = total_clusters;
			}
			
			int temp_array[cluster_size[parent]]; //Array used for linear sorting of the cluster assign
			for(i=0;i<cluster_size[parent];i++)
				temp_array[i]=-1;
			for(i=0;i<cluster_size[parent];i++){
				if(cluster_assign[cluster_start[parent]+i]==parent)
					temp_array[i] = 0;
				else
					temp_array[i] = 1;
			}
			
			int x=0, y=cluster_size[parent]-1,temp1=0;
			double temp2=0.0;
			
			//<---- Sorting the cluster assign along with data points parallely ---->
			while(x<y){
				while(temp_array[x]==0)
					x++;
				while(temp_array[y]==1)
					y--;
				if(x<y){
					temp1 = temp_array[x];
					temp_array[x] = temp_array[y];
					temp_array[y] = temp1;
					temp1 = cluster_assign[cluster_start[parent]+x];
					cluster_assign[cluster_start[parent]+x] = cluster_assign[cluster_start[parent]+y];
					cluster_assign[cluster_start[parent]+y] = temp1;
					for(i=0;i<dim;i++){
						temp2 = data[cluster_start[parent]+x*dim+i];
						data[cluster_start[parent]+x*dim+i] = data[cluster_start[parent]+y*dim+i];
						data[cluster_start[parent]+y*dim+i] = temp2;
					}
					x++;
					y--;
				}
			}
			
			
			brk_cond = 0; // condition used for breaking infinite loop
			
			//<---- Comparing the current and previous cluster assigns
			for(i=0;i<cluster_size[parent];i++){
				if(temp_assign[i]!=cluster_assign[cluster_start[parent]+i] || temp_assign[i]==-1)
					brk_cond = 1;
			}
			
			for(i=0;i<cluster_size[parent];i++) //Re-assigning the current cluster assign to temp assign used for the next iteration of while loop
				temp_assign[i] = cluster_assign[cluster_start[parent]+i];
			
			//<---- Calculating the start and size of both the clusters ------>
			cluster_size[total_clusters] = cluster_size[parent] - x;
			cluster_size[parent] = x;
			cluster_start[total_clusters] = cluster_start[parent] + cluster_size[parent];
			
			if(brk_cond==0) // Checking the breaking condition
				break;
			
			int empty_cluster = 0, other_cluster = 0;
			
			//<---- Checking for the empty cluster ------>
			if(cluster_size[parent]==0 || cluster_size[total_clusters]==0){
				if(cluster_size[parent]==0){
					empty_cluster = parent;
					other_cluster = total_clusters;
				}
				else{
					empty_cluster = total_clusters;
					other_cluster = parent;
				}
				int far_index1 = -1;
				double max_dist1 = LONG_MIN;
				
				//<---- Finding the farthest point from cluster centroid in the non-empty cluster ----->
				for(i=0;i<cluster_size[other_cluster];i++){
					double temp_dist = 0.0;
					for(j=0;j<dim;j++)
						temp_dist = temp_dist + pow(data[cluster_start[other_cluster]+i*dim+j]-cluster_centroid[other_cluster][j],2);
					temp_dist = sqrt(temp_dist);
					if(temp_dist>max_dist1){
						max_dist1 = temp_dist;
						far_index1 = cluster_start[other_cluster]+i;
					}
				}
				
				//<----- Assigning the farthest point to the empty cluster centroid ------>
				for(i=0;i<dim;i++)
					cluster_centroid[empty_cluster][i] = data[far_index1+i];
			}
			
			//<-------- Loop executes if there is no empty cluster ------>
			else{
				
				//<------- Re-calculating both cluster centroids ------>
				for(i=0;i<dim;i++){
					double temp_sum=0.0;
					for(j=0;j<cluster_size[parent];j++)
						temp_sum = temp_sum + data[cluster_start[parent]+j*dim+i];
					temp_sum = temp_sum / cluster_size[parent];
					cluster_centroid[parent][i] = temp_sum;
				}
				for(i=0;i<dim;i++){
					double temp_sum=0.0;
					for(j=0;j<cluster_size[total_clusters];j++)
						temp_sum = temp_sum + data[cluster_start[total_clusters]+j*dim+i];
					temp_sum = temp_sum / cluster_size[total_clusters];
					cluster_centroid[total_clusters][i] = temp_sum;
				}
			}
		}
		
		//<----- Calculating SSE for both the clusters ------->
		double max_sse_dist = LONG_MIN;
		int sse_cluster = -1;
		for(i=0;i<total_clusters;i++){
			double total_dist = 0.0;
			for(j=0;j<cluster_size[i];j++){
				double temp_dist = 0.0;
				for(x=0;x<dim;x++)
					temp_dist = temp_dist + pow(cluster_centroid[i][x]-data[cluster_start[i]+j*dim+x],2);
				total_dist = total_dist + temp_dist;
			}
			if(total_dist>max_sse_dist){
				max_sse_dist = total_dist;
				sse_cluster = i;
			}
		}
		
		parent = sse_cluster; //Assigning the cluster with highest SSE to parent cluster
		total_clusters++;
	}
	kmeans(dim, ndata, data, k, cluster_size, cluster_start, cluster_radius, cluster_centroid, cluster_assign);
}

void kmeans(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double *cluster_radius, double **cluster_centroid, int *cluster_assign){
	int i=0, j=0, x=0, brk_cond=0, temp_assign[ndata];
	double *query; // Stores the query point
	query = (double*) malloc(dim * sizeof(double));
	for(i=0;i<k;i++)
		cluster_size[i]=0; //Initializing all the cluster sizes to zero
	for(i=0;i<ndata;i++)
		temp_assign[i] = cluster_assign[i]; //Initializing the temp assign to the initial cluster assign
	
	//<-------- Infinite loop executes until the current and previous cluster assigns are equal -------->
	while(1){
		
		//<---------- Assigning all the points to clusters by using the distances between point and cluster centroids
		for(i=0;i<ndata;i++){
			double min_dist = LONG_MAX;
			int min_cluster=-1;
			for(j=0;j<k;j++){
				double temp_dist = 0.0;
				for(x=0;x<dim;x++)
					temp_dist = temp_dist + pow(data[i*dim+x]-cluster_centroid[j][x],2);
				temp_dist = sqrt(temp_dist);
				if(temp_dist<min_dist){
					min_dist = temp_dist;
					min_cluster = j;
				}
			}
			cluster_assign[i] = min_cluster;
			cluster_size[min_cluster]++;
		}
		
		brk_cond=0; //Breaking condition to break the infinite loop
		
		//<-------- Comparing the current and previous cluster assigns -------->
		for(i=0;i<ndata;i++){
			if(temp_assign[i]!=cluster_assign[i])
				brk_cond=1;
		}
		
		if(brk_cond==0) //Checking for the breaking condition
			break;
		
		//<------- Re-calculating all the cluster starts
		cluster_start[0] = 0;
		for(i=1;i<k;i++)
			cluster_start[i] = cluster_start[i-1] + cluster_size[i-1];
		
		int cluster_empty = -1;
		
		//<------- Checking for the empty cluster ------->
		for(i=0;i<k;i++){
			if(cluster_size[i]==0)
				cluster_empty = i;
			break;
		}
		
		int empty_index=-1;
		
		//<------ Loop executes if empty cluster is found -------->
		if(cluster_empty!=-1){
			double max_dist = LONG_MIN;
			
			//<------ Finding the minimum of maximum distances from each point to each cluster centroid ------->
			for(i=0;i<ndata;i++){
				double min_dist = LONG_MAX;
				for(j=0;j<k;j++){
					if(j==cluster_empty)
						continue;
					double temp_dist = 0.0;
					for(x=0;x<dim;x++)
						temp_dist = temp_dist + pow((data[i*dim+x]-cluster_centroid[j][x]),2);
					temp_dist = sqrt(temp_dist);
					if(temp_dist<min_dist)
						min_dist = temp_dist;
				}
				if(min_dist>max_dist){
					max_dist = min_dist;
					empty_index = i;
				}
			}
			
			//<------ Assigning the point to empty cluster centroid -------->
			for(i=0;i<dim;i++)
				cluster_centroid[cluster_empty][i] = data[empty_index+i];
		}
		
		//<--------- Loop executes if empty cluster is not found -------->
		else{
			int counter = 0, temp_assign = -1;
			double temp_data = -1.0;
			
			//<------- Swapping the cluster assigns and corresponding data points -------->
			for(i=0;i<k;i++){
				for(j=counter;j<ndata;j++){
					if(cluster_assign[j]==i){
						temp_assign = cluster_assign[counter];
						cluster_assign[counter] = cluster_assign[j];
						cluster_assign[j] = temp_assign;
						for(x=0;x<dim;x++){
							temp_data = data[counter*dim+x];
							data[counter*dim+x] = data[j*dim+x];
							data[j*dim+x] = temp_data;
						}
					}
					counter++;
				}
			}
		}
		
		//<------ Assigning the current cluster centroids to temp assign used in next iteration ------>
		for(i=0;i<ndata;i++)
			temp_assign[i] = cluster_assign[i];
		
		//<----- Initializing all the cluster centroids to zero -------->
		for(i=0;i<k;i++)
				for(j=0;j<dim;j++)
					cluster_centroid[i][j]=0.0;
		
		//<-------- Re-calculating all the cluster centroids --------->
		for(i=0;i<k;i++){
			for(j=0;j<dim;j++){
				double temp_sum = 0.0;
				for(x=0;x<cluster_size[i];x++)
					temp_sum = temp_sum + data[cluster_start[i]+x*dim+j];
				cluster_centroid[i][j] = temp_sum / cluster_size[i];
			}
		}
	}
	
	//<--------- Calculating the cluster radius ------->
	for(i=0;i<k;i++){
		double max_dist = LONG_MIN;
		for(j=0;j<cluster_size[i];j++){
			double temp_dist = 0.0;
			for(x=0;x<dim;x++)
				temp_dist = temp_dist + pow((cluster_centroid[i][x]-data[cluster_start[i]+j*dim+x]),2);
			temp_dist = sqrt(temp_dist);
			if(temp_dist > max_dist)
				max_dist = temp_dist;
		}
		cluster_radius[i] = max_dist;
	}
	
	//<------- Calling the search query for ten random points -------->
	for(i=0;i<10;i++){
		for(j=0;j<dim;j++)
			query[j] = rand() % 101;
		search_kmeans(dim, ndata, data, k, cluster_size, cluster_start, cluster_radius, cluster_centroid, query);
	}
}
		
void search_kmeans(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double *cluster_radius, double **cluster_centroid, double *query){
	int i=0, j=0, visits=0, min_cluster=-1, cluster_visited[k];
	double min_dist=LONG_MAX, cluster_dist[k];
	for(i=0;i<k;i++){
		cluster_dist[i]=-1.0;
		cluster_visited[i]=0;
	}
	
	//<----- Calculating the cluster distances from point to cluster boundary ------>
	for(i=0;i<k;i++){
		double cal_dist=0.0;
		for(j=0;j<dim;j++)
			cal_dist = cal_dist + pow((query[j]-cluster_centroid[i][j]),2);
		cal_dist = sqrt(cal_dist);
		cal_dist = abs(cal_dist - cluster_radius[i]);
		cluster_dist[i] = cal_dist;
	}
	
	//<----- Loop executes until min cluster distance if found ------->
	while(1){
		int brk_cond=0;
		double min=LONG_MAX;
		
		//<----- Finding the min cluster distance ------>
		for(i=0;i<k;i++){
			if(cluster_visited[i]==1)
				continue;
			else if(cluster_visited[i]==0){
				if(cluster_dist[i]<min_dist){
					if(cluster_dist[i]<min){
						min=cluster_dist[i];
						min_cluster = i;
					}
					brk_cond++;
				}
			}
		}
		
		if(brk_cond==0) // Breaking condition for the while loop
			break;
		
		//<-------- Searching for the min distance inside the cluster having least cluster distance ------->
		for(i=0;i<cluster_size[min_cluster];i++){
			double cal_dist=0.0;
			for(j=0;j<dim;j++)
				cal_dist = cal_dist + pow((query[j]-data[cluster_start[min_cluster]+i*dim+j]),2);
			cal_dist = sqrt(cal_dist);
			if(cal_dist<min_dist || min_dist==-1.0)
				min_dist = cal_dist;
			visits++;
		}
		cluster_visited[min_cluster]=1;
	}
	
	printf("Visits = %d\n",visits);
	printf("Minimum Distance = %lf\n",min_dist);
}