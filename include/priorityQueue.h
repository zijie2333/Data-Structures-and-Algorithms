#ifndef PRIORITYQUEUE_H_INCLUDED
#define PRIORITYQUEUE_H_INCLUDED

#include <iostream>

using namespace std;
class box;

template<class Type>
class priorityQueue{
    friend class box;
public:
    priorityQueue(int capacity=100){
        array = new Type[capacity]; maxSize = capacity;currentSize=0;
    }
    priorityQueue(const Type data[],int size);
    ~priorityQueue(){ delete []array;}

    bool isEmpty()const{return currentSize==0;}
    void enQueue(const Type&x);
    Type deQueue();
    Type getHeap()const {return array[1];}
    int findMin(const Type&x)const;
    void traverse()const;
    void decreaseKey(int i,const Type&value);
private:
    int currentSize;
    Type*array;
    int maxSize;

    void doubleSpace();
    void buildHeap();
    void precolateDown(int hole);
    void swap(int x1,int x2){ Type tmp = array[x1];array[x1]=array[x2];array[x2]=tmp;}
};

template<class Type>
void priorityQueue<Type>::enQueue(const Type&x)
{
    if(currentSize == maxSize-1)doubleSpace();
    int hole = ++currentSize;
    for(;hole>1&&x<array[hole/2];hole/=2)
        array[hole] = array[hole/2];
    array[hole]=x;
}
template<class Type>
Type priorityQueue<Type>::deQueue()
{
    Type minItem = array[1];
    array[1] = array[currentSize--];
    precolateDown(1);
    return minItem;
}
template<class Type>
void priorityQueue<Type>::precolateDown(int hole)
{
    int child;
    Type tmp = array[hole];

    for(;hole*2<=currentSize;hole=child)
    {
        child = hole*2;
        if(child!=currentSize && array[child+1]<array[child])child++;
        if(array[child]<tmp)array[hole]=array[child];
        else break;
    }
    array[hole] = tmp;
}
template<class Type>
void priorityQueue<Type>::buildHeap()
{
    for(int i=currentSize/2;i>0;i--)
    {
        precolateDown(i);
    }
}
template<class Type>
priorityQueue<Type>::priorityQueue(const Type*items,int size):maxSize(size+10),currentSize(size)
{
    array = new Type[maxSize];
    for(int i=0;i<size;i++)
        array[i+1]=items[i];
    buildHeap();
}
template<class Type>
void priorityQueue<Type>::doubleSpace()
{
    Type*tmp =array;
    maxSize *=2;
    array = new Type[maxSize];
    for(int i=0;i<=currentSize;++i)array[i]=tmp[i];
    delete []tmp;
}
template<class Type>
int priorityQueue<Type>::findMin(const Type&x)const
{
    Type min;int minID;
    bool first=true;
    for(int i=0;i<currentSize;i++)
        if(array[i]>x)
        {
            if(first){min = array[i];minID=i;first=false;}
            else if(array[i]<min){ min = array[i];minID=i;}
        }
    return minID;
}
template<class Type>
void priorityQueue<Type>::traverse()const
{
    for(int i=1;i<currentSize;i++)
        cout << array[i]<<' ';
}
template<class Type>
void priorityQueue<Type>::decreaseKey(int i,const Type&value)
{
    if(i>currentSize || i<=0)return;
    array[i]-=value;
    // 一定可以取代target的位置
    while(array[i]<array[i/2])
    {
        int target = findMin(array[i]);
        swap(i,target);
    }
}
#endif // PRIORITYQUEUE_H_INCLUDED
