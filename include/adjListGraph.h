#ifndef ADJLISTGRAPH_H_INCLUDED
#define ADJLISTGRAPH_H_INCLUDED

#include <queue>
#include <stack>
#include "DisjointSet.h"
#include "priorityQueue.h"
#include "string"

using namespace std;

template<class TypeOfVer, class TypeOfEdge>
class adjListGraph {
public:
	adjListGraph(int vSize, const TypeOfVer d[]);
	bool insert(int u, int v, TypeOfEdge w);  // u,v顶点间插入权值w的边
	bool insert2(int u, int v, TypeOfEdge w); // 双向边
	bool remove(int u, int v);               // 删除u,v顶点之间的边
	bool exist(int u, int v)const;
	void dfs()const;
	void bfs()const;
	void EulerCircuit(TypeOfVer start);
	vector<vector<TypeOfVer> > topSort()const;  // 拓扑排序
	void find(const TypeOfVer&start, int M);

	void kruskal()const;
	void mykruskal(TypeOfVer x)const;
	void prim(TypeOfEdge noEdge)const;
	void myprim(TypeOfEdge noEdge)const;

	void unweightedShortDistance(TypeOfVer start, TypeOfEdge noEdge)const;
	void Dijkstra(TypeOfVer start,TypeOfEdge noEdge)const;
	void weightNegative(TypeOfVer start, TypeOfEdge noEdge)const;
	void DijkstraLeastNodes(TypeOfVer start, TypeOfEdge noEdge)const;
	void DijkstraPriorityQueue(TypeOfVer start, TypeOfEdge noEdge)const;

	~adjListGraph();

private:
	// 边结点
	struct edgeNode {
		int end;
		TypeOfEdge weight;
		edgeNode* next;
		edgeNode() { next = NULL; }
		edgeNode(int e, TypeOfEdge w, edgeNode*n = NULL)
		{
			end = e; weight = w; next = n;
		}
		bool operator<(const edgeNode&rp)const { return weight<rp.weight; }
	};
	// 顶点结点
	struct verNode {
		TypeOfVer ver;
		edgeNode*head;
		verNode(edgeNode*h = NULL) { head = h; }
	};
	// 欧拉回路结点
	struct EulerNode {
		int NodeNum;
		EulerNode *next;
		EulerNode(int v) { NodeNum = v; next = NULL; }
	};
	// 深度-路径结点
	struct pathNode {
		int nodeNum;
		//int depth;
		pathNode *next;
		pathNode(int v, pathNode*_next = NULL) :nodeNum(v), next(_next) {}
	};
	// 边类型
	struct edge {
		int beg, end;
		TypeOfEdge w;
		bool operator<(const edge&rp)const { return w<rp.w; }
		edge() {}
		edge(int _beg, int _end, const TypeOfEdge &_w) :beg(_beg), end(_end), w(_w) {}
	};
	struct nodeDistance{
		int node;
		int distance;
		nodeDistance(){}
		nodeDistance(int n,int dis):node(n),distance(dis){}
		bool operator<(const nodeDistance& nd)const { return distance < nd.distance; }
	};
	void dfs(int start, bool visited[])const;
	verNode* clone()const;  //复制图
	EulerNode* EulerCircuit(int start, EulerNode*&end);
	void find(int start, int len, int *visitState, pathNode* curPath);
	void showPath(pathNode* path)const;

	void printPath(int start, int end, int prev[])const;

	verNode* verList;
	int Vers;
	int Edges;


};

template<class TypeOfVer, class TypeOfEdge>
adjListGraph< TypeOfVer, TypeOfEdge>::adjListGraph(int vSize, const TypeOfVer d[])
{
	Vers = vSize; Edges = 0;
	verList = new verNode[vSize];
	for (int i = 0; i<Vers; ++i)verList[i].ver = d[i];
}
template<class TypeOfVer, class TypeOfEdge>
adjListGraph< TypeOfVer, TypeOfEdge>::~adjListGraph()
{
	edgeNode*p;
	for (int i = 0; i<Vers; ++i)
	{
		while ((p = verList[i].head) != NULL)
		{
			verList[i].head = p->next;
			delete p;
		}
	}
	delete[]verList;
}
template<class TypeOfVer, class TypeOfEdge>
bool adjListGraph< TypeOfVer, TypeOfEdge>::insert(int u, int v, TypeOfEdge w)
{
	verList[u].head = new edgeNode(v, w, verList[u].head);
	++Edges;
	return true;
}
template<class TypeOfVer, class TypeOfEdge>
bool adjListGraph< TypeOfVer, TypeOfEdge>::insert2(int u, int v, TypeOfEdge w) // 双向边
{
	return insert(u, v, w) && insert(v, u, w);
}
template<class TypeOfVer, class TypeOfEdge>
bool adjListGraph< TypeOfVer, TypeOfEdge>::remove(int u, int v)
{
	edgeNode*tmp = verList[u].head;
	edgeNode*q;
	if (tmp == NULL)return false;
	if (tmp->end == v)
	{
		verList[u].head = tmp->next;
		delete tmp; --Edges;
		return true;
	}
	while (tmp->next != NULL && tmp->next->end != v)tmp = tmp->next;
	if (tmp->next == NULL)return false;
	q = tmp->next; tmp->next = q->next; delete q;
	--Edges;

	return true;
}

template<class TypeOfVer, class TypeOfEdge>
bool adjListGraph< TypeOfVer, TypeOfEdge>::exist(int u, int v)const
{
	edgeNode*p = verList[u].head;
	for (; p->next != NULL; p = p->next)
		if (p->end == v)return true;
	return false;
}

template<class TypeOfVer, class TypeOfEdge>
void adjListGraph< TypeOfVer, TypeOfEdge>::dfs()const
{
	bool *visited = new bool[Vers];
	for (int i = 0; i<Vers; i++)visited[i] = false;
	cout << "当前图像的深度优先遍历序列:" << endl;
	for (int i = 0; i<Vers; i++)
	{
		if (visited[i])continue;
		dfs(i, visited);
		cout << endl;
	}
}
template<class TypeOfVer, class TypeOfEdge>
void adjListGraph< TypeOfVer, TypeOfEdge>::dfs(int start, bool visited[])const
{
	edgeNode*p = verList[start].head;
	cout << verList[start].ver << "\t";
	visited[start] = true;
	while (p != NULL) {
		if (!visited[p->end])dfs(p->end, visited);
		p = p->next;
	}
}
template<class TypeOfVer, class TypeOfEdge>
void adjListGraph< TypeOfVer, TypeOfEdge>::bfs()const
{
	bool *visited = new bool[Vers];
	int curNode;
	queue<int>q;
	edgeNode* edge;
	for (int i = 0; i<Vers; i++)visited[i] = false;

	cout << "当前图的广度优先遍历序列为:" << endl;

	for (int i = 0; i<Vers; i++) {
		if (visited[i] == true)continue;
		q.push(i);
		while (!q.empty())
		{
			curNode = q.front(); q.pop();
			if (visited[curNode])continue;
			cout << verList[curNode].ver << '\t';
			visited[curNode] = true;
			edge = verList[curNode].head;
			while (edge != NULL)
			{
				if (visited[edge->end] == false)q.push(edge->end);
				edge = edge->next;
			}

		}
		cout << endl;
	}
}
template<class TypeOfVer, class TypeOfEdge>
typename adjListGraph<TypeOfVer, TypeOfEdge>::verNode* adjListGraph< TypeOfVer, TypeOfEdge>::clone()const
{
	verNode* tmp = new verNode[Vers];
	edgeNode*p;
	for (int i = 0; i<Vers; i++)
	{
		tmp[i].ver = verList[i].ver;
		p = verList[i].head;
		while (p != NULL)
		{
			tmp[i].head = new edgeNode(p->end, p->weight, tmp[i].head);
			p = p->next;
		}
	}
	return tmp;
}
template<class TypeOfVer, class TypeOfEdge>
typename adjListGraph<TypeOfVer, TypeOfEdge>::EulerNode*
adjListGraph< TypeOfVer, TypeOfEdge>::EulerCircuit(int start, typename adjListGraph<TypeOfVer, TypeOfEdge>::EulerNode*&end)
{
	EulerNode*beg;
	int nextNode;

	beg = end = new EulerNode(start);
	while (verList[start].head != NULL)
	{
		nextNode = verList[start].head->end;
		remove(start, nextNode);
		remove(nextNode, start);
		start = nextNode;
		end->next = new EulerNode(start);
		end = end->next;
	}
	return beg;
}
template<class TypeOfVer, class TypeOfEdge>
void adjListGraph< TypeOfVer, TypeOfEdge>::EulerCircuit(TypeOfVer start)
{
	EulerNode* beg, *end, *p, *q, *tb, *te;
	int numOfDegree;
	edgeNode*r;
	verNode*tmp;

	if (Edges == 0) { cout << "不存在欧拉回路" << endl; return; }
	for (int i = 0; i<Vers; i++)
	{
		numOfDegree = 0;
		r = verList[i].head;
		while (r != NULL) { numOfDegree++; r = r->next; }
		if (numOfDegree == 0 || numOfDegree % 2) { cout << "不存在欧拉回路" << endl; return; }
	}

	int i = 0;
	for (; i<Vers; i++)
		if (verList[i].ver == start)break;
	if (i == Vers) { cout << "起始节点不存在" << endl; return; }

	tmp = clone();

	beg = EulerCircuit(i, end);
	while (true)
	{
		p = beg;
		while (p->next != NULL)
			if (verList[p->next->NodeNum].head != NULL)break;
			else p = p->next;
			if (p->next == NULL)break; // 所有边都访问过了
			q = p->next;
			tb = EulerCircuit(q->NodeNum, te);   // 找到新回路
			te->next = q->next;                 // 插入原路径
			p->next = tb;
			delete q;  // 新路径中包含了q 故将原来的q释放
	}
	delete[]verList;
	verList = tmp;

	cout << "欧拉回路是:" << endl;
	while (beg != NULL)
	{
		cout << verList[beg->NodeNum].ver << '\t';
		p = beg; beg = beg->next;
		delete p;
	}
	cout << endl;
}
template<class TypeOfVer, class TypeOfEdge>
vector<vector<TypeOfVer> > adjListGraph< TypeOfVer, TypeOfEdge>::topSort()const
{
	queue<int> q;
	int *degree = new int[Vers];
	for (int i = 0; i<Vers; i++)degree[i] = 0;
	vector<vector<TypeOfVer> >semester;

	edgeNode *tmp;
	for (int i = 0; i < Vers; i++) //计算入度
	{
		tmp = verList[i].head;
		while (tmp != NULL)
		{
			degree[tmp->end] += tmp->weight;
			tmp = tmp->next;
		}
	}
	for (int i = 0; i < Vers; i++)if (degree[i] == 0)q.push(i);
	q.push(-1);
	semester.push_back(vector<TypeOfVer>());
	int sem = 0;
	while (!q.empty())
	{
		int cur = q.front(); q.pop();
		if (cur == -1) {
			// 这个学期能学的课都学完了 下个学期可以学的课都入队了 增加标识
			sem++;
			if (!q.empty()) { q.push(-1); semester.push_back(vector<TypeOfVer>()); }
			continue;
		}
		semester[sem].push_back(cur);
		for (edgeNode*p = verList[cur].head; p != NULL; p = p->next) {
			degree[p->end] -= p->weight;
			if (degree[p->end] == 0)q.push(p->end);
		}

	}

	return semester;
}
class ERRORCLASS {};
#define GOOD 0
#define VISITED 1
#define FULL 2

template<class TypeOfVer, class TypeOfEdge>
void adjListGraph< TypeOfVer, TypeOfEdge>::find(const TypeOfVer&start, int M)
{
	int i = 0;
	for (; i < Vers; i++)
		if (verList[i].ver == start)break;
	if (i == Vers)throw ERRORCLASS();
	int *visitState = new int[Vers];
	for (int j = 0; j < Vers; j++)visitState[j] = GOOD;
	pathNode* beg = NULL;
	find(i, M, visitState, beg);
}
template<class TypeOfVer, class TypeOfEdge>
void adjListGraph< TypeOfVer, TypeOfEdge>::find(int start, int len, int *visitState, pathNode* curPath) //能不能修改curPath?
{
	visitState[start] |= VISITED;
	pathNode *tmpPath = new pathNode(start, curPath);

	if (len == 1)
	{
		for (edgeNode*p = verList[start].head; p != NULL; p = p->next)
		{
			if (!(visitState[p->end] & VISITED))// 未访问过
			{
				pathNode* resPath = new pathNode(p->end, tmpPath);
				showPath(resPath);  // 输出path
				visitState[p->end] |= FULL;
				delete resPath;
			}
		}
	}
	else { // 检查所有邻接点
		edgeNode* p = verList[start].head;
		for (; p != NULL; p = p->next)
		{
			if (visitState[p->end] & VISITED)continue;
			if (!(visitState[p->end] & FULL))
				find(p->end, len - 1, visitState, tmpPath);
		}
		if (p == NULL)//所有邻接点都 访问过或者满 (必然)
			visitState[start] = FULL; // 设为 未访问 满
	}
	delete tmpPath;
}
template<class TypeOfVer, class TypeOfEdge>
void adjListGraph<TypeOfVer, TypeOfEdge>::showPath(pathNode* path)const
{
	stack<int> s;
	for (; path != NULL; path = path->next)
		s.push(path->nodeNum);
	int cur;
	while (!s.empty())
	{
		cur = s.top(); s.pop();
		cout << verList[cur].ver;
		if (!s.empty()) cout << " --> ";
	}
	cout << endl;
}
template<class TypeOfVer, class TypeOfEdge>
void adjListGraph<TypeOfVer, TypeOfEdge>::kruskal()const
{
	int edgeAccepted = 0;
	DisjointSet ds(Vers);
	priorityQueue<edge> q;

	for (int i = 0; i < Vers; i++)
		for (edgeNode*p = verList[i].head; p != NULL; p = p->next)
			if (i < p->end)q.enQueue(edge(i, p->end, p->weight));
	while (edgeAccepted < Vers - 1)
	{
		edge cur = q.deQueue();
		if (ds.Find(cur.beg) != ds.Find(cur.end))
		{
			edgeAccepted++;
			ds.Union(ds.Find(cur.beg), ds.Find(cur.end));
			cout << '(' << verList[cur.beg].ver << ',' << verList[cur.end].ver << ')' << '\t';
		}
	}
	cout << endl;
}
template<class TypeOfVer, class TypeOfEdge>
void adjListGraph<TypeOfVer, TypeOfEdge>::mykruskal(TypeOfVer x)const
{
	edgeNode *edges = new edgeNode[Vers];
	bool *flag = new bool[Vers];

	int edgeAccepted = 0;
	DisjointSet ds(Vers);
	priorityQueue<edge> q;

	int start = 0;
	for (; start < Vers; start++)if (verList[start].ver == x)break;
	if (start == Vers)return;

	for (int i = 0; i < Vers; i++)
	{
		flag[i] = false;
		for (edgeNode*p = verList[i].head; p != NULL; p = p->next)
			if (i < p->end)q.enQueue(edge(i, p->end, p->weight));
	}
	while (edgeAccepted < Vers - 1)
	{
		edge cur = q.deQueue();
		if (ds.Find(cur.beg) != ds.Find(cur.end))
		{
			edgeAccepted++;
			ds.Union(ds.Find(cur.beg), ds.Find(cur.end));
			//cout << '(' << verList[cur.beg].ver << ',' << verList[cur.end].ver << ')' << endl;
			edges[cur.beg].next = new edgeNode(cur.end, cur.w, edges[cur.beg].next);
			edges[cur.end].next = new edgeNode(cur.beg, cur.w, edges[cur.end].next);
		}
	}
	// get minimum tree from edges
	stack<int>s;
	edgeNode*p, *tmp;

	s.push(start);
	flag[start] = true;
	while (!s.empty())
	{
		start = s.top();
		for (p = edges[start].next; p != NULL; p = p->next)
		{
			if (!flag[p->end])
			{
				flag[p->end] = true;
				s.push(p->end);
				cout << '(' << verList[start].ver << ',' << verList[p->end].ver << ')' << '\t';
				break; // 找到一个可访问的邻结点 直接入栈 结束
			}
		}
		if (p == NULL)s.pop(); //这个结点访问完毕 出栈
	}

	for (int i = 0; i<Vers; i++)
		for (p = edges[i].next; p != NULL;)
		{
			tmp = p;
			p = p->next;
			delete tmp;
		}
	delete[]flag;
	delete[]edges;
}
template<class TypeOfVer, class TypeOfEdge>
void adjListGraph<TypeOfVer, TypeOfEdge>::prim(TypeOfEdge noEdge)const
{
	bool *flag = new bool[Vers];
	TypeOfEdge* lowCost = new TypeOfEdge[Vers];
	int *startNode = new int[Vers];

	edgeNode*p;
	TypeOfEdge min;
	int start, i, j;
	for (i = 0; i < Vers; i++)
	{
		flag[i] = false;
		lowCost[i] = noEdge;
	}
	start = 0;
	for (i = 1; i < Vers; i++)
	{
		for (p = verList[start].head; p != NULL; p = p->next)
		{
			if (!flag[p->end] && lowCost[p->end] > p->weight) //更新距离
			{
				lowCost[p->end] = p->weight;
				startNode[p->end] = start;
			}
		}
		flag[start] = true;
		min = noEdge;
		for (j = 0; j < Vers; j++)
			if (lowCost[j] < min)
			{
				min = lowCost[j];
				start = j;
			}
		cout << '(' << verList[startNode[start]].ver << ',' << verList[start].ver << ')' << '\t';
		lowCost[start] = noEdge;

	}
	delete[]flag;
	delete[]startNode;
	delete[]lowCost;
}
template<class TypeOfVer, class TypeOfEdge>
void adjListGraph<TypeOfVer, TypeOfEdge>::myprim(TypeOfEdge noEdge)const
{
	bool *flag = new bool[Vers];
	priorityQueue<edge> q;
	edge curEdge;
	int i, start;

	for (i = 0; i < Vers; i++)flag[i] = false;
	int cnt = 1;
	curEdge.end = 0;
	start = 0;
	while (cnt < Vers)
	{
		flag[start] = true;
		for (edgeNode*p = verList[start].head; p != NULL; p = p->next)
			if (!flag[p->end])q.enQueue(edge(start, p->end, p->weight));
		while (flag[curEdge.end]) curEdge = q.deQueue();
		start = curEdge.end;
		cnt++;
		cout << '(' << verList[curEdge.beg].ver << ',' << verList[curEdge.end].ver << ')' << '\t';
	}
}
template<class TypeOfVer, class TypeOfEdge>
void adjListGraph<TypeOfVer, TypeOfEdge>::unweightedShortDistance(TypeOfVer start, TypeOfEdge noEdge)const
{
	queue<int> q;
	TypeOfEdge *distance = new TypeOfEdge[Vers];
	int *prev = new int[Vers];
	int beg, cur;
	edgeNode *tmp;

	for (int i = 0; i<Vers; i++) distance[i] = noEdge;
	for (beg = 0; beg<Vers; beg++)
		if (verList[beg].ver == start)break;
	if (beg == Vers) { cout << start << " is not in the graph!" << endl; return; }

	distance[beg] = 0;
	prev[beg] = beg;
	q.push(beg);

	while (!q.empty())
	{
		cur = q.front(); q.pop();
		for (tmp = verList[cur].head; tmp != NULL; tmp = tmp->next)
			if (distance[tmp->end] == noEdge)
			{
				distance[tmp->end] = distance[cur] + 1;
				prev[tmp->end] = cur;
				q.push(tmp->end);
			}
	}
	for (int i = 0; i<Vers; i++)
	{
		cout << "从" << start << "到" << verList[i].ver << "的路径:" << endl;
		printPath(beg, i, prev);
		cout << endl;
	}
}
template<class TypeOfVer, class TypeOfEdge>
void adjListGraph<TypeOfVer, TypeOfEdge>::printPath(int start, int end, int prev[])const
{
	if (start == end) { cout << verList[start].ver; return; }
	printPath(start, prev[end], prev);
	cout << " -> " << verList[end].ver;
}

template<class TypeOfVer, class TypeOfEdge>
void adjListGraph<TypeOfVer, TypeOfEdge>::Dijkstra(TypeOfVer start, TypeOfEdge noEdge)const
{
	bool *known = new bool[Vers];
	int *prev = new int[Vers];
	TypeOfEdge *distance = new TypeOfEdge[Vers];

	TypeOfEdge min;
	edgeNode *p;

	int sNo, u;
	for (int i = 0; i < Vers; i++) {
		known[i] = false;
		distance[i] = noEdge;
	}
	for (sNo = 0; sNo < Vers; sNo++)
		if (verList[sNo].ver == start)break;
	if (sNo == Vers) { cout << start << " is not in graph!" << endl; return; }

	distance[sNo] = 0;
	prev[sNo] = sNo;
	for (int i = 1; i < Vers; i++)
	{
		min = noEdge;
		for(int j=0;j<Vers;j++)
			if (!known[j] && distance[j] < min)
			{
				min = distance[j];
				u = j;
			}
		known[u] = true;
		for (p = verList[u].head; p != NULL; p = p->next)
		{
			if (!known[p->end] && distance[p->end] > min + p->weight)
			//if ( distance[p->end] > min + p->weight)
			{
				distance[p->end] = min + p->weight;
				prev[p->end] = u;
			}
		}
	}

	for (int i = 0; i < Vers; i++)
	{
		cout << "从" << start << "到" << verList[i].ver << "的路径:" << endl;
		printPath(sNo, i, prev);
		cout << "长度为:" << distance[i] << endl << endl;
	}
}
template<class TypeOfVer, class TypeOfEdge>
void adjListGraph<TypeOfVer, TypeOfEdge>::weightNegative(TypeOfVer start,TypeOfEdge noEdge)const
{
	int *prev = new int[Vers];
	TypeOfEdge *distance = new TypeOfEdge[Vers];
	queue<int> q;
	edgeNode*p;
	int cur;
	int sNo = 0;
	for (int i = 0; i < Vers; i++) distance[i] = noEdge;

	for (; sNo < Vers; sNo++)
		if (verList[sNo].ver == start)break;
	if (sNo == Vers) { cout << start << " is not in graph!" << endl; return; }

	q.push(sNo);
	prev[sNo] = sNo;
	distance[sNo] = 0;
	while (!q.empty())
	{
		cur = q.front(); q.pop();
		for (p = verList[cur].head; p != NULL; p = p->next)
		{
			if (p->weight + distance[cur] < distance[p->end])
			{
				distance[p->end] = p->weight + distance[cur];
				prev[p->end] = cur;
				q.push(p->end);
			}
		}
	}
	for (int i = 0; i < Vers; i++)
	{
		cout << "从" << start << "到" << verList[i].ver << "的路径:" << endl;
		printPath(sNo, i, prev);
		cout << "长度为:" << distance[i] << endl << endl;
	}
}
template<class TypeOfVer, class TypeOfEdge>
void adjListGraph<TypeOfVer, TypeOfEdge>::DijkstraLeastNodes(TypeOfVer start, TypeOfEdge noEdge)const
{
	bool *known = new bool[Vers];
	int *prev = new int[Vers];
	int *nodes = new int[Vers];
	TypeOfEdge *distance = new TypeOfEdge[Vers];

	TypeOfEdge min;
	edgeNode *p;

	int sNo, u;
	for (int i = 0; i < Vers; i++) {
		known[i] = false;
		distance[i] = noEdge;
		nodes[i] = -1;
	}
	for (sNo = 0; sNo < Vers; sNo++)
		if (verList[sNo].ver == start)break;
	if (sNo == Vers) { cout << start << " is not in graph!" << endl; return; }

	distance[sNo] = 0;
	prev[sNo] = sNo;
	nodes[sNo] = 1;
	for (int i = 1; i < Vers; i++)
	{
		min = noEdge;
		for (int j = 0; j<Vers; j++)
			if (!known[j] && distance[j] < min)
			{
				min = distance[j];
				u = j;
			}
		known[u] = true;
		for (p = verList[u].head; p != NULL; p = p->next)
		{
			if (!known[p->end] && distance[p->end] > min + p->weight)
			{
				distance[p->end] = min + p->weight;
				prev[p->end] = u;
				nodes[p->end] = nodes[u] + 1;
			}
			else if (!known[p->end] && distance[p->end] == min + p->weight)
			{
				if (nodes[p->end] > nodes[u] + 1)
				{
					distance[p->end] = min + p->weight;
					prev[p->end] = u;
					nodes[p->end] = nodes[u] + 1;
				}
			}
		}
	}

	for (int i = 0; i < Vers; i++)
	{
		cout << "从" << start << "到" << verList[i].ver << "的路径:" << endl;
		printPath(sNo, i, prev);
		cout << "长度为:" << distance[i] << endl << endl;
	}
}
template<class TypeOfVer, class TypeOfEdge>
void adjListGraph<TypeOfVer, TypeOfEdge>::DijkstraPriorityQueue(TypeOfVer start, TypeOfEdge noEdge)const
{
	bool *known = new bool[Vers];
	int *prev = new int[Vers];
	TypeOfEdge *distance = new TypeOfEdge[Vers];

	edgeNode *p;

	int sNo;
	for (int i = 0; i < Vers; i++) {
		known[i] = false;
		distance[i] = noEdge;
	}
	for (sNo = 0; sNo < Vers; sNo++)
		if (verList[sNo].ver == start)break;
	if (sNo == Vers) { cout << start << " is not in graph!" << endl; return; }

	distance[sNo] = 0;
	prev[sNo] = sNo;

	nodeDistance cur;
	priorityQueue<nodeDistance> pq;
	pq.enQueue(nodeDistance(sNo,distance[sNo]));
	for (int i = 1; i < Vers; i++)
	{
		do {
			cur = pq.deQueue();
		} while (known[cur.node]);
		known[cur.node] = true;
		for (p = verList[cur.node].head; p != NULL; p = p->next)
		{
			if (!known[p->end] && distance[p->end] > cur.distance + p->weight)
			{
				distance[p->end] = cur.distance + p->weight;
				prev[p->end] = cur.node;
				pq.enQueue(nodeDistance(p->end, distance[p->end]));
			}
		}
	}

	for (int i = 0; i < Vers; i++)
	{
		cout << "从" << start << "到" << verList[i].ver << "的路径:" << endl;
		printPath(sNo, i, prev);
		cout << "长度为:" << distance[i] << endl << endl;
	}
}
#endif // ADJLISTGRAPH_H_INCLUDED
