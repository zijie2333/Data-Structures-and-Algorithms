#ifndef DISJOINTSET_H_INCLUDED
#define DISJOINTSET_H_INCLUDED

class DisjointSet{
private:
    int size;
    int*parent;
public:
    DisjointSet(int s);
    ~DisjointSet(){delete[]parent;}
    void Union(int root1,int root2);
    int Find(int x);

};


#endif // DISJOINTSET_H_INCLUDED
