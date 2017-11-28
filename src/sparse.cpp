#include "sparse.h"
#include <cassert>
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
  Ap=new int [cap];
  Ai=new int [2*cap];
  TAp=new int[cap];
  TAi=new int[2*cap];
  Ax=new double[2*cap];
  TAx=new double[2*cap];
    
  y=new double[cap];
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
  delete index;
    
  Ap=new int [cap];
  Ai=new int [2*cap];
  TAp=new int[cap];
  TAi=new int[2*cap];
  Ax=new double[2*cap];
  TAx=new double[2*cap];
  
  y=new double[cap];
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
  
  vector<SparseMatrixElem> ttmps;
  
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
        ttmps.push_back(te);
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
    ttmps.push_back(te);
  }
  dim = max(r, c)+1;
  for (int j = c; j < dim; j++) {
    Ap[j + 1] = nz;
  }

  
  stable_sort(ttmps.begin(), ttmps.end());
  c = 0;
  r = 0;
  nz = 0;
  
  for (size_t i = 0; i < ttmps.size(); i++) {
    int tempColumn = ttmps[i].column;
    r = ttmps[i].row;

    if(c!=tempColumn){
      for (int j = c; j < tempColumn; j++) {
        TAp[j + 1] = nz;
      }
      c=tempColumn;
    }
    TAx[nz]=ttmps[i].value;
    TAi[nz++] = r;

  }

  for (int j = c; j < dim; j++) {
    TAp[j + 1] = nz;
  }
  
  nonzero=nz;
  
  computableConnectedComponent();
}

SparseSolver::~SparseSolver() {

  if(NULL!=NULL){
    delete[] Ap;
    delete[] Ai;
    delete[] TAp;
    delete[] TAi;
    delete[] Ax;
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
  I.clear();
  J.clear();
  fill(bb.begin(), bb.begin()+dim, false);


  for(int i=0; i< dim; i++){
    row[i]=TAp[i+1]-TAp[i];
  }
  
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
  bool state=true;
  
  while(state){
    state=false;
    
    vector<int>::iterator it=I.begin();
    while(it!=I.end()){
      if((1==row[*it])&& (!bb[*it]) ){
        int start=TAp[*it];
        int end=TAp[*it+1];
        int c=-1;
        for(int i=start; i< end; i++){
          int p=TAi[i];
          if(!deleteColumn[p]){
            c=p;
          }
        }
        
        state=true;
        deleteColumn[c]=true;
        
        start=Ap[c];
        end=Ap[c+1];
        
        for(int i=start; i< end; i++){
          int r=Ai[i];
          row[r]-=1;
        }

        it=I.erase(it);

      }else{
        it++;
      }
    }
  }
  vector<int> TJ(J);
  J.clear();
  for(vector<int>::iterator it=TJ.begin(); it!= TJ.end(); it++){
    if(!deleteColumn[*it]){
      J.push_back(*it);
    }
  }

#ifdef DEBUG
  cout<<"final size: "<<I.size()<<endl;
#endif
  
  sort(I.begin(), I.end());
  sort(J.begin(), J.end());
}

void SparseSolver::tminComputableProjection(const SparseVector& b,
                                            vector<int>& I,
                                            vector<int>& J)  {
  I.clear();
  J.clear();
  fill(bb.begin(), bb.begin()+dim, false);



  for(int i=0; i< dim; i++){
    row[i]=Ap[i+1]-Ap[i];
  }
    
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
  bool state=true;
  
  while(state){
    state=false;
    
    vector<int>::iterator it=I.begin();
    while(it!=I.end()){
      if((1==row[*it])&& (!bb[*it]) ){
        int start=Ap[*it];
        int end=Ap[*it+1];
        int c=-1;
        for(int i=start; i< end; i++){
          int p=Ai[i];
          if(!deleteColumn[p]){
            c=p;
          }
        }
        
        state=true;
        deleteColumn[c]=true;
        
        start=TAp[c];
        end=TAp[c+1];
        
        for(int i=start; i< end; i++){
          int r=TAi[i];
          row[r]-=1;
        }

        it=I.erase(it);

      }else{
        it++;
      }
    }
  }

  vector<int> TJ(J);
  J.clear();
  for(vector<int>::iterator it=TJ.begin(); it!= TJ.end(); it++){
    if(!deleteColumn[*it]){
      J.push_back(*it);
    }
  }

#ifdef DEBUG
  cout<<"final size: "<<I.size()<<endl;
#endif
  
  sort(I.begin(), I.end());
  sort(J.begin(), J.end());
}
bool SparseSolver::locSolver(SparseVector& sB, int* ap, int* ai, double* ax,
                             vector<int>& vecI, vector<int>& vecJ)  {
  int smallest = vecI.front();
  int bigest = vecI.back();

  cs* B = cs_spalloc(0, 0, 1, 1, 1);

  for (size_t j = 0; j < vecJ.size(); j++) {
    int k = vecJ[j];
    int start = ap[k];
    int end = ap[k + 1];
    for (int l = start; l < end; l++) {
      int p = ai[l];
      if (p >= smallest && p <= bigest) {
        int i = binfind(vecI.begin(), vecI.end(), p);
        if (i > -1) {
          cs_entry(B, i, j, ax[l]);
        }
      }
    }
  }


  fill(y, y + vecJ.size(), 0);
  for (size_t j = 0; j < sB.locs.size(); j++) {
    int k = sB.locs[j];
    int i = binfind(vecI.begin(), vecI.end(), k);
    assert(i > -1);
    y[i] = sB.values[j];
  }

  cs* M = cs_compress(B);

  cs_spfree(B);

  klu_symbolic* Symbolic;
  klu_numeric* Numeric;
  klu_common Common;
  klu_defaults(&Common);
  int n = M->n;

  for (int i = 0; i <= n; i++) {
    p[i] = M->p[i];
  }
  int num = M->p[n];

  for (int i = 0; i < num; i++) {
    index[i] = M->i[i];
  }

  Symbolic = klu_analyze(M->n, p, index, &Common);
  Numeric = klu_factor(p, index, M->x, Symbolic, &Common);

  int re = klu_solve(Symbolic, Numeric, M->n, 1, y, &Common);

  klu_free_symbolic(&Symbolic, &Common);
  klu_free_numeric(&Numeric, &Common);

  cs_spfree(M);
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
  
  for (csi i = 0; i < dim; i++) {

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

  for (csi i = 0; i < dim; i++) {

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
