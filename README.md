#ICP 3D registration with ANN

##introduction
This lab is aim to realize the 3D point cloud registration between two dense sets of points. Through Iterative Closest Points (ICP) Algorithm, we can estimate transformation parameters. The points in {Set 1} using estimated parameters to close {Set 2}. 

In the ICP algorithm, we need to search the nearest neighbor of each point. If without acceleration, it calculates the distance for each point in {Set 1} and {Set 2}. Time complexity is O(nm), n = number of {Set 1}, m = number of {Set 2}. So we have to find ways to reduce the running time. 

ANN is a library to exact and approximate nearest neighbor searching. We can use ANN to build tree data structure. After that, nearest neighbor searching is will be easier. If with ANN, time complexity is O(N). 

This lab shows two version codes. One is the basic version, one is the ANN acceleration version. And it shows their error, running time and average CPU. Through these, we can compare their complexity.

##file structure
Include ann library.

main function: icp_registration.cpp

You can find detials in the report: TP_ICP_Report_ZixinZOU.docx