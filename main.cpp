
#include <iostream>
#include <string>
#include "./include/adjListGraph.h"

using namespace std;

int main()
{
	cout << "��14�µ�6�� Dijkstra�����ȼ����иĽ�" << endl;
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
	g.insert(3, 5, 8);  // ԭͼ
	//g.insert(3, 5, -8); // ��Ȩֵʹ��
	//g.insert(3, 5, 5);   // ��С���Dijkstraʹ��
	g.insert(3, 6, 4);
	g.insert(3, 4, 2);
	g.insert(4, 6, 6);
	g.insert(6, 5, 1);
	//cout << endl << "unweightedShortDistance:" << endl;
	//g.unweightedShortDistance(string("V2"), 10000);
	//cout << endl << "Dijkstra:" << endl;
	//g.Dijkstra(string("V1"), 10000);
	//cout << endl << "weightNegative(V2��ʼ):" << endl;
	//g.weightNegative(string("V2"), 10000);
	//cout << endl << "DijkstraLeastNodes(��V1��ʼ):" << endl;
	//g.DijkstraLeastNodes(string("V1"), 10000);
	//cout << "��V1��V5������·�� �ֱ�Ϊ��\nV1 -> V3 -> V5 ����8" << endl;
	//cout << "V1 -> V3 -> V6 ->V5 ����8" << endl;
	//cout << "V1 -> V4 -> V6 ->V5 ����17" << endl;
	//cout << "DijkstraLeastNodes ѡ���˵�һ��·����������С�ҽ��������" << endl << endl;


	cout << endl << "Dijkstra:" << endl;
	g.Dijkstra(string("V1"), 10000);
	cout << endl << "DijkstraPriorityQueue:" << endl;
	g.DijkstraPriorityQueue(string("V1"), 10000);
	return 0;
}
