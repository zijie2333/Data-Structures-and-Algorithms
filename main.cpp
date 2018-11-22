
#include <iostream>
#include <string>
#include "./include/adjListGraph.h"

using namespace std;

int main()
{
	cout << "第14章第6题 Dijkstra的优先级队列改进" << endl;
	//int test[] = { 1,2,3,4,5,6 };
	string test[] = { "V0","V1","V2","V3","V4","V5","V6" };
	adjListGraph<string, int> g(7, test);

	g.insert(0, 1, 2);
	g.insert(0, 3, 1);
	g.insert(1, 4, 10);
	g.insert(1, 3, 3);
	g.insert(2, 0, 4);
	g.insert(2, 5, 5);
	g.insert(3, 2, 2);
	g.insert(3, 5, 8);  // 原图
	//g.insert(3, 5, -8); // 负权值使用
	//g.insert(3, 5, 5);   // 最小结点Dijkstra使用
	g.insert(3, 6, 4);
	g.insert(3, 4, 2);
	g.insert(4, 6, 6);
	g.insert(6, 5, 1);
	//cout << endl << "unweightedShortDistance:" << endl;
	//g.unweightedShortDistance(string("V2"), 10000);
	//cout << endl << "Dijkstra:" << endl;
	//g.Dijkstra(string("V1"), 10000);
	//cout << endl << "weightNegative(V2开始):" << endl;
	//g.weightNegative(string("V2"), 10000);
	//cout << endl << "DijkstraLeastNodes(从V1开始):" << endl;
	//g.DijkstraLeastNodes(string("V1"), 10000);
	//cout << "从V1到V5有三条路径 分别为：\nV1 -> V3 -> V5 长度8" << endl;
	//cout << "V1 -> V3 -> V6 ->V5 长度8" << endl;
	//cout << "V1 -> V4 -> V6 ->V5 长度17" << endl;
	//cout << "DijkstraLeastNodes 选择了第一条路径，长度最小且结点数最少" << endl << endl;


	cout << endl << "Dijkstra:" << endl;
	g.Dijkstra(string("V1"), 10000);
	cout << endl << "DijkstraPriorityQueue:" << endl;
	g.DijkstraPriorityQueue(string("V1"), 10000);
	return 0;
}
