#include "sparse.h"
#include <cassert>
#include<deque>
#include<iostream>
#include "graph.h"
#include "graphalg.hpp"
#include "klu.h"
#include "util.h"
namespace raptor {
namespace sparse {
using namespace std;

subMatrix subMatrix::operator+(const subMatrix& other) const {
  subMatrix re = other;
  re.add(other);

  return re;
}

void subMatrix::add(const subMatrix& other) {
  rows.insert(rows.end(), rows.begin(), rows.end());
  columns.insert(columns.end(), columns.begin(), columns.end());
  sort(rows.begin(), rows.end());
  sort(columns.begin(), columns.end());

#ifdef DEBUG
  for (size_t i = 0; i + 1 < rows.size(); i++) {
    assert(rows[i] != rows[i + 1]);
  }
  for (size_t i = 0; i + 1 < columns.size(); i++) {
    assert(columns[i] != columns[i + 1]);
  }
#endif
}

SparseSolver::SparseSolver() {
  nonzero=0;
  y=NULL;
  p=index=NULL;
  cap=100;
  dim=0;
  Ap=new int [cap+1];
  Ai=new int [2*cap];
  TAp=new int[cap+1];
  TAi=new int[2*cap];
  Ax=new double[2*cap];
  TAx=new double[2*cap];
    
  y=new double[2*cap];
  p=new int[cap+1];
  index= new int [2*cap];
  bb.resize(cap);
  row.resize(cap);
  deleteColumn.resize(cap);

}
void SparseSolver::reScale(int s){

  cap=s;
  cap*=1.3;

  delete[] Ap;
  delete[] Ai;
  delete[] TAp;
  delete[] TAi;
  delete[] Ax;
  delete[] TAx;
  
  delete[] y;
  delete[] p;
  delete[] index;
    
  Ap=new int [cap+1];
  Ai=new int [2*cap];
  TAp=new int[cap+1];
  TAi=new int[2*cap];
  Ax=new double[2*cap];
  TAx=new double[2*cap];
  
  y=new double[2*cap];
  p=new int[cap+1];
  index= new int [2*cap];
    
  bb.resize(cap);
  row.resize(cap);
  deleteColumn.resize(cap);

  
}
void SparseSolver::update( vector<SparseMatrixElem>& elements) {

  if((int)elements.size()>2*cap || dim+1>cap ){
    int s=max((int)elements.size()/2, dim+1)+1;
    reScale(s);
  }
  
  vector<SparseMatrixElem> columnElems;
  
  stable_sort(elements.begin(), elements.end());
  int c = 0;
  int r = 0;
  double value = 0;
  Ap[0]=TAp[0]=0;
  int nz = 0;
  for (size_t i = 0; i < elements.size(); i++) {
    int tempColumn = elements[i].column;
    int tempRow = elements[i].row;
    if (c == tempColumn && r == tempRow) {
      value += elements[i].value;
    } else {
      if (fabs(value) > 1e-6) {
        Ai[nz] = r;
        Ax[nz] = value;
        nz++;

        SparseMatrixElem te;
        te.row=c;
        te.column=r;
        te.value=value;
        columnElems.push_back(te);
      }

      if (c != tempColumn) {
        for (int j = c; j < tempColumn; j++) {
          Ap[j + 1] = nz;
        }
      }
      c = tempColumn;
      r = tempRow;
      value = elements[i].value;
    }
  }
  if (fabs(value) > 1e-6) {
    Ai[nz] = r;
    Ax[nz] = value;
    nz++;

    SparseMatrixElem te;
    te.row=c;
    te.column=r;
    te.value=value;
    columnElems.push_back(te);
  }
  dim = max(r, c)+1;
  for (int j = c; j < dim; j++) {
    Ap[j + 1] = nz;
  }

  
  sort(columnElems.begin(), columnElems.end());
  c = 0;
  r = 0;
  nz = 0;
  
  for (size_t i = 0; i < columnElems.size(); i++) {
    int tempColumn = columnElems[i].column;
    r = columnElems[i].row;

    if(c!=tempColumn){
      for (int j = c; j < tempColumn; j++) {
        TAp[j + 1] = nz;
      }
      c=tempColumn;
    }
    TAx[nz]=columnElems[i].value;
    TAi[nz++] = r;

  }

  for (int j = c; j < dim; j++) {
    TAp[j + 1] = nz;
  }
  
  nonzero=nz;
  
  computableConnectedComponent();
}

SparseSolver::~SparseSolver() {

  if(NULL!=Ap){
    delete[] Ap;
    delete[] Ai;
    delete[] TAp;
    delete[] TAi;
    delete[] Ax;
    delete[] TAx;
  }
  
  if(NULL!=y){
    delete[] y;
  }
  if(NULL!=p){
    delete[] p;
  }
  if(NULL!=index){
    delete[] index;
  }
  
}

void SparseSolver::minComputableProjection(const SparseVector& b,
                                           vector<int>& I,
                                           vector<int>& J)  {

  fill(bb.begin(), bb.begin()+dim, false);

  
  set<int> passSubs;
  
  for (vector<int>::const_iterator it = b.locs.begin(); it != b.locs.end();
       it++) {

    passSubs.insert(subIndex[*it]);
    bb[*it]=true;

  }

  for (set<int>::iterator it = passSubs.begin(); it != passSubs.end(); it++) {
    I.insert(I.end(), subs[*it].rows.begin(), subs[*it].rows.end());
    J.insert(J.end(), subs[*it].columns.begin(), subs[*it].columns.end());
  }
  
#ifdef DEBUG
  cout<<"orignal size: "<<I.size()<<endl;
#endif


  fill(deleteColumn.begin(), deleteColumn.begin()+dim, false);

  vector<int> waitS;
  for(vector<int>::iterator it=I.begin(); it!= I.end(); it++){
    row[*it]=TAp[*it+1]-TAp[*it];
    if((1==row[*it])&& (!bb[*it]) ){
      waitS.push_back(*it);
    }
  }
  vector<int> secondS;
  

  while(!waitS.empty()){
    for(vector<int>::iterator it=waitS.begin(); it!=waitS.end(); it++){
    
      int k=*it;

      int start=TAp[k];
      int end=TAp[k+1];
      int c=-1;
      for(int i=start; i< end; i++){
        int p=TAi[i];
        if(!deleteColumn[p]){
          c=p;
        }
      }
      deleteColumn[c]=true;
        
      start=Ap[c];
      end=Ap[c+1];
        
      for(int i=start; i< end; i++){
        int r=Ai[i];
        row[r]-=1;
        if((1==row[r]) &&(!bb[r]) ){
          secondS.push_back(r);
        }
      }
    }
    waitS.swap(secondS);
    secondS.clear();
  }

  vector<int> TJ;
  TJ.swap(J);

  for(vector<int>::iterator it=TJ.begin(); it!= TJ.end(); it++){
    if(!deleteColumn[*it]){
      J.push_back(*it);
    }
  }
  vector<int> TI;
  TI.swap(I);

  for(vector<int>::iterator it=TI.begin(); it!= TI.end(); it++){
    if(row[*it]>0){
      I.push_back(*it);
    }
  }
  

#ifdef DEBUG
  cout<<"final size: "<<I.size()<<endl;
#endif
  

}

void SparseSolver::tminComputableProjection(const SparseVector& b,
                                            vector<int>& I,
                                            vector<int>& J)  {

  fill(bb.begin(), bb.begin()+dim, false);

    
  set<int> passSubs;
  for (vector<int>::const_iterator it = b.locs.begin(); it != b.locs.end();
       it++) {
    passSubs.insert(subIndex[*it + dim]);
    bb[*it]=true;
  }

  for (set<int>::iterator it = passSubs.begin(); it != passSubs.end(); it++) {
    I.insert(I.end(), subs[*it].columns.begin(), subs[*it].columns.end());
    J.insert(J.end(), subs[*it].rows.begin(), subs[*it].rows.end());
  }

#ifdef DEBUG
  cout<<"orignal size: "<<I.size()<<endl;
#endif

  
  fill(deleteColumn.begin(), deleteColumn.begin()+dim, false);

  vector<int> waitS;
  for(vector<int>::iterator it=I.begin(); it!= I.end(); it++){
    row[*it]=Ap[*it+1]-Ap[*it];
    if((1==row[*it])&& (!bb[*it]) ){
      waitS.push_back(*it);
    }
  }
  vector<int> secondS;
  while(!waitS.empty()){
    for(vector<int>::iterator it=waitS.begin(); it!=waitS.end(); it++){
      int k=*it;

      int start=Ap[k];
      int end=Ap[k+1];
      int c=-1;
      for(int i=start; i< end; i++){
        int p=Ai[i];
        if(!deleteColumn[p]){
          c=p;
        }
      }

      deleteColumn[c]=true;
        
      start=TAp[c];
      end=TAp[c+1];
        
      for(int i=start; i< end; i++){
        int r=TAi[i];
        row[r]-=1;
        if((1==row[r])&& (!bb[r]) ){
          secondS.push_back(r);
        } 
      }
    }
    waitS.swap(secondS);
    secondS.clear();
  }

  vector<int> TJ;
  TJ.swap(J);
  for(vector<int>::iterator it=TJ.begin(); it!= TJ.end(); it++){
    if(!deleteColumn[*it]){
      J.push_back(*it);
    }
  }

  vector<int> TI;
  TI.swap(I);
  for(vector<int>::iterator it=TI.begin(); it!= TI.end(); it++){
    if(row[*it]>0){
      I.push_back(*it);
    }
  }

#ifdef DEBUG
  cout<<"final size: "<<I.size()<<endl;
#endif
  

}
bool SparseSolver::locSolver(SparseVector& sB, int* ap, int* ai, double* ax,
                             vector<int>& vecI, vector<int>& vecJ)  {
  
  fill(row.begin(), row.begin()+dim, -1);
  for(size_t i=0; i< vecI.size(); i++){
    row[vecI[i]]=i;
  }
  int n=vecI.size();
  int nz=0;
  p[0]=0;
  for(int j=0; j< n; j++){
    int k = vecJ[j];
    int start = ap[k];
    int end = ap[k + 1];
    for (int l = start; l < end; l++) {
      int p = ai[l];
        
      int i =row[p];
      if (i > -1) {
        index[nz]=i;
        y[nz]=ax[l];
        nz++;
      }
    }
    p[j+1]=nz;
  }

  klu_symbolic* Symbolic;
  klu_numeric* Numeric;
  klu_common Common;
  klu_defaults(&Common);

  Symbolic = klu_analyze(n, p, index, &Common);
  Numeric = klu_factor(p, index, y, Symbolic, &Common);


  fill(y, y + vecJ.size(), 0);
  for (size_t j = 0; j < sB.locs.size(); j++) {
    int k = sB.locs[j];
    int i = row[k];
    assert(i > -1);
    y[i] = sB.values[j];
  }


  int re = klu_solve(Symbolic, Numeric, n, 1, y, &Common);

  klu_free_symbolic(&Symbolic, &Common);
  klu_free_numeric(&Numeric, &Common);


  if (1 != re) {
    return false;
  }
  /**
   * lifting
   *
   */

  sB.clear();
  for (size_t i = 0; i < vecJ.size(); i++) {
    if (fabs(y[i]) > 1e-6) {
      int k = vecJ[i];
      sB.locs.push_back(k);
      sB.values.push_back(y[i]);
    }
  }
  return true;
}

void SparseSolver::computableConnectedComponent() {
  if (dim>0) {
    subs.clear();
    
    subIndex.resize(2 * dim);
    fill(subIndex.begin(), subIndex.end(), 0);


    fill(bb.begin(), bb.begin()+dim, true);
    fill(deleteColumn.begin(), deleteColumn.begin()+dim, true);
    vector<int> tempR, tempC;
    for (int i = 0; i < dim; i++) {
      if ((bb[i]) && (TAp[i+1]>TAp[i]) ) {
        bb[i]=false;
        subMatrix temp;
        int k=subs.size();
        tempR.push_back(i);
        while(!tempR.empty()){
          temp.rows.insert( temp.rows.end(), tempR.begin(), tempR.end());
         
          for(vector<int>::iterator it=tempR.begin(); it!= tempR.end(); it++){
            int r=*it;
            subIndex[r]=k;
            int start=TAp[r];
            int end=TAp[r+1];
            for(int j=start;  j< end; j++){
              int c=TAi[j];
              if(deleteColumn[c]){
                deleteColumn[c]=false;
                tempC.push_back(c);
              }
            }            
          }

          temp.columns.insert(temp.columns.end(), tempC.begin(), tempC.end());
          tempR.clear();

          for(vector<int>::iterator it=tempC.begin(); it!=tempC.end(); it++){
            int c=*it;
            subIndex[c+dim]=k;
            int start=Ap[c];
            int end=Ap[c+1];
            for(int j=start;  j< end; j++){
              int r=Ai[j];
              if(bb[r]){
                bb[r]=false;
                tempR.push_back(r);
              }
            }            
            
          }
          tempC.clear();
        }

        subs.push_back(temp);
      }
    }
  }
}

bool SparseSolver::locSolver(SparseVector& sB)  {
  if(sB.locs.empty()){
    return true;
  }


  vector<int> vecI, vecJ;

  minComputableProjection(sB, vecI, vecJ);

  assert(vecI.size() == vecJ.size());

#ifdef DEBUG
  cout<<"dimentsion: "<<dim<<" submatrix dimension: "<<vecI.size()<<endl;
#endif

  return locSolver(sB, Ap, Ai, Ax, vecI, vecJ);
}
bool SparseSolver::locSolver(double* b)  {
  SparseVector sB;
  for (int i = 0; i < dim; i++) {
    if (fabs(b[i]) > 1e-8) {
      sB.locs.push_back(i);
      sB.values.push_back(b[i]);
    }
  }
#ifdef DEBUG
  cout<<"sparse vector size: "<<sB.locs.size()<<endl;
#endif
  if (!locSolver(sB)) {
    return false;
  }
  fill(b, b + dim, 0);
#ifdef DEBUG
   cout<<"sparse solution vector size: "<<sB.locs.size()<<endl;
#endif
  for (size_t i = 0; i < sB.locs.size(); i++) {
    b[sB.locs[i]] = sB.values[i];
  }
  return true;
}

bool SparseSolver::tlocSolver(SparseVector& sB)  {
  if(sB.locs.empty()){
    return true;
  }

  vector<int> vecI, vecJ;

  tminComputableProjection(sB, vecI, vecJ);

  assert(vecI.size() == vecJ.size());
  return locSolver(sB, TAp, TAi, TAx, vecI, vecJ);
}

bool SparseSolver::tlocSolver(double* b)  {
  SparseVector sB;
  for (int i = 0; i < dim; i++) {
    if (fabs(b[i]) > 1e-6) {
      sB.locs.push_back(i);
      sB.values.push_back(b[i]);
    }
  }
  if (!tlocSolver(sB)) {
    return false;
  }
  fill(b, b + dim, 0);
  for (size_t i = 0; i < sB.locs.size(); i++) {
    b[sB.locs[i]] = sB.values[i];
  }
  return true;
}

bool SparseSolver::incSolver(const double *initSolution, double *b) {
  
  for (int i = 0; i < dim; i++) {

    int start = TAp[i];
    int end = TAp[i + 1];
    for (int k = start; k < end; k++) {
      int j = TAi[k];
      b[i]-=TAx[k]*initSolution[j];
    }
  }

  if(!locSolver(b)){
    return false;
  }
  for(int i=0; i< dim; i++){
    b[i]+=initSolution[i];
  }
  return true;
  
}


bool SparseSolver::tincSolver(const double *initSolution, double *b) {

  for (int i = 0; i < dim; i++) {

    int start = Ap[i];
    int end = Ap[i + 1];
    for (int k = start; k < end; k++) {
      int j = Ai[k];
      b[i]-=Ax[k]*initSolution[j];
    }
  }

  if(!tlocSolver(b)){
    return false;
  }
  for(int i=0; i< dim; i++){
    b[i]+=initSolution[i];
  }
  return true;
  
}


}
}
