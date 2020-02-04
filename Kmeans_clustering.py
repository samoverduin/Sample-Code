#!/usr/bin/env python3
"""
Author: Sam Overduin
Implementation of the k-means clustering algorithm

"""

# import statements
import random
import math
import statistics

sys_random = random.SystemRandom()

def csv_parser(file):
    """Return list of point coordinates as [[x1,y1,z1,...],[x2,y2,z2,...]]
    
    file: str of file relative location. Expected format:
        The file has a single header line. 
        Each line contains the coordinates for one data point, starting
        with a label. A data point can be specified in arbitrary dimensions.

    Output: List of lists with the coordinates of the points.
    This function does not capture the labels of the data points. In case
    they are needed later, this function should be adjusted. 
    """
    data_points = []
    f = open(file, "r")
    for line in f:
        items = line.strip().split(",")
        try: #will fail on header line in file
            data_points.append(list(map(float, items[1:]))) # first item is the
                                                            # label
        except ValueError: #must be the header
            continue
    f.close()
    return data_points


def calc_euclid_dist(point1, point2):
    """Calculate the Euclidean distance between point 1 & 2

    :param point1: lst or tuple of [dimensions]
    :param point2: lst or tuple of [dimensions]
    :return: euclidean distance
    """
    assert len(point1) == len(point2), "input points not same # of dimensions"
    squared_sum = 0
    for i in range(len(point1)):
        squared_sum += (point1[i]-point2[i])**2
    return math.sqrt(squared_sum)


def init_centroids(point_lst, k):
    """ Returns k random centroids from point_lst

    :param point_lst: list of point tuples [label, dimensions]
    :param k: int of k
    :return: list of k centroids
    """
    centroids = []
    while len(centroids) < k:
        new_point = random.choice(point_lst)
        if new_point not in centroids:
            centroids.append(new_point)
    return centroids


def associate_cluster(point_lst, centroids):
    """find closest centroids to points in point_lst

    :param point_lst:
    :param centroids:
    :return: lst of closest centroids in point_lst given by index [1...k+1]
    """
    cluster_lst = []

    for point in point_lst:
        centroid_dist_lst = []
        for centroid in centroids:
            centroid_dist_lst.append(calc_euclid_dist(point, centroid))
        # append minimum centroid index + 1 (cluster #)
        cluster_lst.append(centroid_dist_lst.index(min(centroid_dist_lst))+1)
    return cluster_lst


def make_centroids(point_lst, cluster_lst):
    """ returns a list of centroids from cluster-assigned points

    :param point_lst: lst of points
    :param cluster_lst: lst of point cluster assignment ints
    :return: lst of centroids per cluster ordered 1,2,3 etc.
    """
    assert len(point_lst) == len(cluster_lst), \
        "point_lst and cluster_lst should be same length"
    centroid_lst = []
    for clust in range(1,len(set(cluster_lst))+1):
        clust_points = [point_lst[i] for i, x in enumerate(cluster_lst)
                        if x == clust]
        centroid = []
        for i in range(len(clust_points[0])):
            average = sum([x[i] for x in clust_points])/len(clust_points)
            centroid.append(average)
        centroid_lst.append(centroid)
    return centroid_lst


def find_clusters(point_lst, k):
    """Returns point list with k clusters assigned

    :param point_lst: list of coordinates in d space
    :param k: int of clusters to make
    :return: -list of ints showing points assigned to cluster [1...k+1]
             -lst of centroids
    """
    assert len({len(x) for x in point_lst}) == 1, \
        "not all points same # of dimensions"
    new_centroids = init_centroids(point_lst, k)
    old_centroids = None
    i = 0
    while new_centroids != old_centroids:
        i += 1
        # in case it keeps flipflopping
        if i > 500:
            x = None
            while x not in ['y', 'n']:
                x = input('500 iterations have passed, do you want to break'
                          ' (y/n)?')
            if x == 'y':
                break
            else:
                i = 0
                continue
        cluster_lst = associate_cluster(point_lst, new_centroids)
        # sorts out loss of clusters
        if len(set(cluster_lst)) != k:
            missing_k_s = list(set(range(1, k+1)) - set(cluster_lst))
            # print('clusters added:', missing_k_s)
            # getting indices in the for loop gave the same result both loops
            # also, random choices gave me the same 2 options at the beginning
            # solved by adding more samples + set
            indices = list(set(random.choices(list(range(len(cluster_lst))),
                                              k=len(missing_k_s)+20)))
            for i, new_clust in enumerate(missing_k_s):
                cluster_lst[indices[i]] = new_clust

        old_centroids = new_centroids
        new_centroids = make_centroids(point_lst, cluster_lst)

    return cluster_lst, new_centroids, i


def calc_wgss(point_lst, cluster_lst, centroids):
    """Calculates total WGSS

    :param point_lst: list of coordinates in d space
    :param cluster_lst: lst of point cluster assignment ints
    :param centroids: lst of centroid coordinates
    :return: total WGSS
    """
    wgss_lst = []
    wgss = 0
    for clust in range(1, len(set(cluster_lst))+1):
        clust_points = [point_lst[i] for i, x in enumerate(cluster_lst)
                        if x == clust]
        centroid = centroids[clust-1]
        for point in clust_points:
            wgss += calc_euclid_dist(point, centroid)**2
        wgss_lst.append(wgss/len(clust_points))
        wgss = 0
    return sum(wgss_lst)


def calc_bgss(point_lst, cluster_lst, centroids):
    """Calculates BGSS

    :param point_lst: list of coordinates in d space
    :param cluster_lst: lst of point cluster assignment ints
    :param centroids: lst of centroid coordinates
    :return: total BGSS
    """
    bgss = 0

    all_point_lst = [1]*len(point_lst)
    center_point = make_centroids(point_lst, all_point_lst)[0]

    i = 1
    for centroid in centroids:
        elements = len([point_lst[j] for j, x in enumerate(cluster_lst)
                        if x == i])
        bgss += elements * calc_euclid_dist(centroid, center_point)**2
        i += 1
    return bgss


def cluster_wrapper(point_lst, k, iterations=1, print_enabled=False,
                    return_all_w_values=False):
    """ Finds clusters and can return best (by W value) cluster_lst and W-value

    :param point_lst: list of coordinates in d space
    :param k: int of clusters to make
    :param iterations: number of clusters to make before selecting best one
    :param print_enabled: bool to enable printing of iterations to convergence
    :param return_all_w_values: bool to return only all W-values

    :return: best cluster_lst by W-value, centroids
    """
    w_val_lst = []
    cluster_lst_lst = []
    convergence_lst = []
    for i in range(iterations):
        cluster_lst, centroids, steps_to_convergence = find_clusters(point_lst, k)
        wgss = calc_wgss(point_lst, cluster_lst, centroids)
        bgss = calc_bgss(point_lst, cluster_lst, centroids)
        convergence_lst.append(steps_to_convergence)
        cluster_lst_lst.append(cluster_lst)
        w_val_lst.append(wgss/bgss)

    best_indice = w_val_lst.index(min(w_val_lst))
    if print_enabled:
        print('Clusterings performed:', iterations)
        print('Iterations till convergence (mean):',
              statistics.mean(convergence_lst))
        print('Best cluster found:')
        cluster_lst = cluster_lst_lst[best_indice]
        for i in range(0, len(cluster_lst), 25):
            print(cluster_lst[i:i + 25])
        print('Mean W-value:', statistics.mean(w_val_lst))

    if return_all_w_values:
        return w_val_lst
    else:
        return cluster_lst_lst[best_indice], w_val_lst[best_indice]


def cluster_resampling(point_lst, k, method='j', iterations=1, print_info=True):
    """Resamples dataset and returns w-values

    :param point_lst: list of coordinates in d space
    :param k: int of clusters to make
    :param method: str of 'b', 'j' or 's' to respectively bootstrap, jackknife
                    or subsample.
    :param iterations: number of clusters to make before selecting best one per
                        new sampled distribution
    :param print_info: option to print mean, stdev etc.
    :return: w-values of samples
    """
    assert method in ['b', 'j', 's'], "method should be in 'b', 'j' or 's' to " \
                                      "respectively bootstrap, jackknife or " \
                                      "subsample."
    w_values = []
    # resample over observations, not genes so transpose
    t_point_lst = list(map(list, zip(*point_lst)))
    if method == 'b':
        for i in range(len(t_point_lst)):
            # sample with replacement
            t_bootstrap = [random.choice(t_point_lst) for x in t_point_lst]
            bootstrap = list(map(list, zip(*t_bootstrap)))  # transpose
            cluster_lst, w_val = cluster_wrapper(bootstrap, k, iterations,
                                                 print_enabled=False)
            w_values.append(w_val)
    elif method == 'j':
        for i in range(len(t_point_lst)):
            # skip 1 observation every time.
            t_jackknife = [x for x in t_point_lst if x != t_point_lst[i]]
            jackknife = list(map(list, zip(*t_jackknife)))  # transpose
            cluster_lst, w_val = cluster_wrapper(jackknife, k, iterations,
                                                 print_enabled=False)
            w_values.append(w_val)
    elif method == 's':
        for i in range(len(t_point_lst)):
            # lists aren't hashable so tuples used. same as bootstrap but made
            # as a set so no replacement
            t_subsample = list(set(tuple(random.choice(t_point_lst))
                                   for x in t_point_lst))
            subsample = list(map(list, zip(*t_subsample)))  # transpose
            cluster_lst, w_val = cluster_wrapper(subsample, k, iterations,
                                                 print_enabled=False)
            w_values.append(w_val)
    if print_info:
        print('Resampling method:', method)
        print('Resamples done:', len(w_values))
        print('Average W-value:', statistics.mean(w_values))
        print('St. dev. W-value:', statistics.stdev(w_values))
    return w_values


def main():
    try:
        large_set_1 = csv_parser("misc/LargeSet_1.csv")
    except IOError:
        print("File not found: misc/LargeSet_1.csv did you download all files from repository?")
        quit()

    # Find best k-value clustering in LargeSet_1.csv
    print("# finding best k-value clustering in LargeSet_1.csv")
    k_vals_to_test = range(2, 7)
    w_val_lst = []
    cluster_lst_lst = []
    for k in k_vals_to_test:
        print('# For k =', k)
        cluster_lst, w_val = cluster_wrapper(large_set_1, k, iterations=1,
                                             print_enabled=False)
        print('Best W-value:', w_val)
        w_val_lst.append(w_val)
        cluster_lst_lst.append(cluster_lst)
        for i in range(1, k + 1):
            print('Cluster', i, 'contains', cluster_lst.count(i), 'elements')

    print('\n# Clustering for LargeSet_1.csv:')
    print('a) Lowest W-value:', min(w_val_lst))
    best_k = k_vals_to_test[w_val_lst.index(min(w_val_lst))]
    print('b) Overal best k-value:', best_k)
    cluster_lst = cluster_lst_lst[w_val_lst.index(min(w_val_lst))]
    print('c) Best cluster:')
    for i in range(0, len(cluster_lst), 25):
        print(cluster_lst[i:i + 25])
    print('d) Stats of 20 iteration around k =', best_k)
    w_values = cluster_wrapper(large_set_1, best_k, iterations=20,
                               return_all_w_values=True)
    print('Average W-value:', statistics.mean(w_values))
    print('St. dev. W-value:', statistics.stdev(w_values))

    # Cluster resampling:
    print('\n# Cluster resampling bootstrap:')
    cluster_resampling(large_set_1, 5, 'b', print_info=True)
    print('\n# Cluster resampling jackknife:')
    cluster_resampling(large_set_1, 5, 'j', print_info=True)
    print('\n# Cluster resampling subsample:')
    cluster_resampling(large_set_1, 5, 's', print_info=True)
    w_values = cluster_wrapper(large_set_1, 5, iterations=24,
                               return_all_w_values=True)
    print('\nWithout resampling:')
    print('Average W-value:', statistics.mean(w_values))
    print('St. dev. W-value:', statistics.stdev(w_values))




if __name__ == "__main__":
    main()
