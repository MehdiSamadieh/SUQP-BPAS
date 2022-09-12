

#include <gmpxx.h>
#include "SUQP_Support.h"
#include "../../Documents/svn/BPAS/Main.new/include/ExpressionTree/ExprTreeNode.hpp"
#include "../../Documents/svn/BPAS/Main.new/include/ExpressionTree/ExpressionTree.hpp"
#include "../../Documents/svn/BPAS/Main.new/include/Symbol/Symbol.hpp"
#include "../../Documents/svn/BPAS/Main.new/tests/MapleTestTool/MapleTestTool.hpp"
#include <math.h>
#include <iostream>
#include <sstream>

AltArr_t* buildRandomPoly(time_t seed, int nterms, unsigned long int coefBound, degree_t sparsity, int includeNeg) {
  if (nterms <= 0) {
    return NULL;
  }

  static int initRand = 0;
  static gmp_randstate_t R_STATE;
  if (!initRand) {
    time_t t = seed;//time(NULL);
    srand(t);

    gmp_randinit_default (R_STATE);
    gmp_randseed_ui(R_STATE, seed);

    fprintf(stderr, "seed: %lu\n", t);
    initRand = 1;
  }

  degree_t maxTotalDeg = sparsity * (nterms);
  degree_t maxUniDeg = (degree_t) ceil(pow(maxTotalDeg, 1.0 ) - 1);
////////////////////////////////////////////////////////
  /*Node* head;
  Node* tail;

  Node* n = (Node*) malloc(sizeof(Node));
  n->next = NULL;*/
  
  AltArr_t* n = (AltArr_t*) malloc(sizeof(AltArr_t));
  
  n->size=0;
  int alloc=nterms;
  n->elems = (AAElem_t*) malloc(sizeof(AAElem_t)*nterms);
//////////////////////////////////////////////////////////////////////////////

  degree_t degs = 0;// (degrees_t) calloc(nvar, sizeof(degree_t));
  if (sparsity == 2) {
    //if dense always include the constant
    degs = 0;
  } else {
    degs = rand() % 2;    //50/50 chance of first term being a constant;
  }
  degree_t lastDegs = degs;

  mpz_t mpzNum, mpzDen;
  mpz_init(mpzNum);
  mpz_init(mpzDen);

  while(mpz_sgn(mpzNum) == 0) {
    mpz_urandomb(mpzNum, R_STATE, coefBound);
  }
  while(mpz_sgn(mpzDen) == 0) {
    mpz_urandomb(mpzDen, R_STATE, coefBound);
  }

  // long int coef_l = rand() % coefBound;
  // if (coef_l == 0) {
    // coef_l = 1l;
  // }

  if (includeNeg && rand() % 2) {
    //50/50 chance of being negative
    mpz_neg(mpzNum, mpzNum);
    // coef_l *= -1;
  }
  n->elems[(n->size)].deg = degs;
  mpq_init(n->elems[n->size].coef);

  mpz_set(mpq_numref(n->elems[n->size].coef), mpzNum);
  mpz_set(mpq_denref(n->elems[n->size].coef), mpzDen);
  mpq_canonicalize(n->elems[n->size].coef);
  // mpq_set_si(n->coef, coef_l, 1ul);
  // n->coef = coef_ul;
 n->size +=1;
  --nterms;

  degree_t step = 0;
  for(int i = 0; i < nterms; ++i) {
    // n = (Node*) malloc(sizeof(Node));

   // n->next = NULL;
    
    //n->elems;
    step = 0;
    while(step == 0) {
      step = rand() % sparsity;
    }
    ////////////
    degs = lastDegs + step;
//////////////////////////////////////
    mpz_urandomb(mpzNum, R_STATE, coefBound);
    while(mpz_sgn(mpzNum) == 0) {
      mpz_urandomb(mpzNum, R_STATE, coefBound);
    }
    mpz_urandomb(mpzDen, R_STATE, coefBound);
    while(mpz_sgn(mpzDen) == 0) {
      mpz_urandomb(mpzDen, R_STATE, coefBound);
    }

    // coef_l = 0l; 
    // while (coef_l == 0) {
      // coef_l = rand() % coefBound;
    // }
    if (includeNeg && rand() % 2) {
      mpz_neg(mpzNum, mpzNum);
      // coef_l *= -1;
    }

    n->elems[n->size].deg = degs;
    lastDegs = degs;
    mpq_init(n->elems[n->size].coef);
    mpz_set(mpq_numref(n->elems[n->size].coef), mpzNum);
    mpz_set(mpq_denref(n->elems[n->size].coef), mpzDen);
    mpq_canonicalize(n->elems[n->size].coef);
    n->size+=1;
    // mpq_set_si(n->coef, coef_l, 1l);
    // n->coef = coef_ul;

   // tail->next = n;
    //tail = n;
  }
 int itera;
  //Now reverse the poly so it is in decreasing order;
  if ((n->size)%2==0)
  {
     itera = ((n->size)-1)/2;

  }
  else
  {
    itera=((n->size)/2)-1;
  }
  for(int i=0;i<=itera;++i)
  { 
     
   degree_t temp_degree= n->elems[i].deg;

    ratNum_t temp_coef;
    mpq_init(temp_coef);
    mpq_set(temp_coef, n->elems[i].coef);
      
   
   n->elems[i].deg= n->elems[(n->size)-1-i].deg;
  

    mpq_init(n->elems[i].coef);
    mpq_set(n->elems[i].coef,n->elems[(n->size)-1-i].coef);

   
   n->elems[(n->size)-1-i].deg= temp_degree;

    mpq_init(n->elems[(n->size)-1-i].coef);
    mpq_set(n->elems[(n->size)-1-i].coef,temp_coef);
    
  }

  /*
  Node* cur = head;
  Node* prev = head;
  Node* next = head->next;
  cur->next = NULL;
  while(next != NULL) {
    cur = next;
    next = cur->next;
    cur->next = prev;
    prev = cur;
  }
  head = cur;*/

  return n;
}
//////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Construct an ExprTreeNode (well, multiple) which represents a single
 * term. That is, a coefficient and a monomial.
 */
ExprTreeNode* exprTreeNodeFromAAElem(const AAElem_t* n, const Symbol& sym) {
    if (n == NULL) {
        return new ExprTreeNode(0l);
    }

    degree_t degs = n->deg;
    ExprTreeNode* t = new ExprTreeNode(mpq_class(n->coef));
    
    
        degree_t deg = n->deg;

        if (deg > 1) {
            //TODO pass symbol directly to expression tree
            ExprTreeNode* var = new ExprTreeNode(sym);
            ExprTreeNode* num = new ExprTreeNode(deg);
            ExprTreeNode* exp = ExprTreeNode::combineExprTreeNodes(var, num, EXPR_EXP);
            t = ExprTreeNode::combineExprTreeNodes(t, exp, EXPR_MULT);
        } else if (deg == 1) {
            ExprTreeNode* var = new ExprTreeNode(sym);
            t = ExprTreeNode::combineExprTreeNodes(t, var, EXPR_MULT);
        }
    

    return t;
}
//////////////////////////////////////////////////////////////////////////////////
ExpressionTree convertToExpressionTree(AltArr_t* a)  {
    if (a->size==0) { //define zero
        ExprTreeNode* r = new ExprTreeNode(0l);
        ExpressionTree t(r);
        return t; 
    }
/////////////////////

    
   
    Symbol sym = "x";
    ExprTreeNode* prev = exprTreeNodeFromAAElem(a->elems, sym);
    for (int i = 1; i < a->size; ++i) {
        ExprTreeNode* thisNode = exprTreeNodeFromAAElem(&(a->elems[i]), sym);
        prev = ExprTreeNode::combineExprTreeNodes(prev, thisNode, EXPR_ADD);
    }

   
    return ExpressionTree(prev);
}

//////////////////////////////////////////////////////////////////////// polytostring

static std::string  polyToString(AltArr_t* n) {
  int n_size= n->size;
  if ((n->elems)== NULL) {
    return "0";
  }

  std::stringstream ss;

  bool first = true;
  bool needsMult = false;
  bool isConst = true;
  mpq_t coef;
   int counter=0;
   mpq_t i;
     mpq_init(i);
     mpq_set_d(i,1);
  while (/*node != NULL*/ counter < n_size ) {
   mpq_init(coef);


    mpq_set(coef,n->elems[counter].coef);

    //coef = ratNum_class(n->elems[counter].coef);

    // coef = node->coef;
    isConst = true;
    if (coef < 0) {
      mpq_neg(coef,coef);
      //coef *= -1; 
      ss << " - ";
    } else if (!first) {
      ss << " + ";
    }

    
     
    if (/*coef != 1*/  (mpq_cmp(coef,i))) {
      ss << coef;
      needsMult = true;
    }/////////////////////////////////////////////////// omit loop and deg
    degree_t degs = n->elems[counter].deg;
    
      if (degs == 0) {
      counter =counter +1;

        continue;
      }
      isConst = false;
      if (needsMult) {
        ss << "*";
      }
      ss << "x";
      if (degs > 1) {
        ss << "^" << degs;
      }
      needsMult = true;
    

    /*node = node->next;*/ 
    counter =counter +1;
    first = false;
    needsMult = false;
  }
  
  if (isConst && coef == i) {
    ss << coef;
  }

  return ss.str();
}
/////////////////////////////////////////////////////////////////////////////////////////
void testAdditionAndMultiplication(time_t seed, int nterms, unsigned long int coefBound, degree_t sparsity, int includeNeg) {
   AltArr_t* a = buildRandomPoly(seed, nterms, coefBound, sparsity, includeNeg);
   AltArr_t* b = buildRandomPoly(seed, nterms, coefBound, sparsity, includeNeg);
  // Node* c = addPolynomials(a,b,nvar);

 
  AltArr_t* ca = addPolynomials_AA(a, b);


  ExpressionTree aTree = convertToExpressionTree(a);
  ExpressionTree bTree = convertToExpressionTree(b);
  ExpressionTree cTree = convertToExpressionTree(ca);

  cTree -= (aTree + bTree);

  MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();
    std::string retErr;
    if(mapleTest->testIfZero(cTree, &retErr) == 0) {
      std::cerr << "ADDITION TEST FAILED: " << std::endl;
    std::cerr << "a:= " << polyToString(a) << ";" << std::endl;
    std::cerr << "b:= " << polyToString(b ) << ";" << std::endl;
    std::cerr << "ca:= " << polyToString(ca) << ";" << std::endl;
    std::cerr << "retErr " <<  retErr << std::endl;
    } else {
      std::cerr << "Addition test passed" << std::endl;
    }

    
    freePolynomial_AA(ca);

  ca = multiplyPolynomials_AA(a, b);

  aTree = convertToExpressionTree(a);
  bTree = convertToExpressionTree(b);
  cTree = convertToExpressionTree(ca);

  cTree -= (aTree * bTree);
  if(mapleTest->testIfZero(cTree, &retErr) == 0) {
      std::cerr << "MULTIPLICATION TEST FAILED: " << std::endl;
    std::cerr << "a:= " << polyToString(a) << ";" << std::endl;// SMQP_CPPSUPPORT.H
    std::cerr << "b:= " << polyToString(b ) << ";" << std::endl;
    std::cerr << "c:= " << polyToString(ca) << ";" << std::endl;
    std::cerr << "retErr " <<  retErr << std::endl;
    }
    std::cerr << "Multiplication test passed" << std::endl;

   }

////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
void testDivision(time_t seed, int nterms, unsigned long int coefBound, degree_t sparsity, int includeNeg) {
   AltArr_t* a = buildRandomPoly(seed, nterms, coefBound, sparsity, includeNeg);
   AltArr_t* b = buildRandomPoly(seed, nterms, coefBound, sparsity, includeNeg);
  // Node* c = addPolynomials(a,b,nvar);

    AltArr_t* ca = NULL;
    AltArr_t* r = NULL;


   dividePolynomials_AA(a, b,&ca,&r);


  ExpressionTree aTree = convertToExpressionTree(a);
  ExpressionTree bTree = convertToExpressionTree(b);
  ExpressionTree cTree = convertToExpressionTree(ca);
  ExpressionTree rTree = convertToExpressionTree(r);
 
  std::cerr << "C tree: " << cTree.toMapleString() << std::endl;
  std::cerr << "R tree: " << rTree.toMapleString() << std::endl;

  //cTree -= (aTree / bTree);

  std::vector<std::string> inputs;
  inputs.push_back(aTree.toMapleString());
  inputs.push_back(bTree.toMapleString());
  inputs.push_back("x");

  // std::cerr << "a:= " << aTree.toMapleString() << ";" << std::endl;
  // std::cerr << "b:= " << bTree.toMapleString() << ";" << std::endl;
  // std::cerr << "c:= " << cTree.toMapleString() << ";" << std::endl;
  // std::cerr << "r:= " << rTree.toMapleString() << ";" << std::endl;


    MapleTestTool* mapleTest = MapleTestTool::getMapleTestTool();
    mapleTest->restartMapleKernel();
    std::string retErr;

    if (mapleTest->testProcReturn("quo", inputs, cTree.toMapleString(), &retErr) == 0) {

      std::cerr << "DIVISION TEST FAILED: " << std::endl;
      std::cerr << "a:= " << polyToString(a) << ";" << std::endl;
      std::cerr << "b:= " << polyToString(b ) << ";" << std::endl;
      std::cerr << "ca:= " << polyToString(ca) << ";" << std::endl;
      std::cerr << "retErr " <<  retErr << std::endl;
      exit(1);
    } 

    if (mapleTest->testProcReturn("rem", inputs, rTree.toMapleString(), &retErr) == 0) {
      std::cerr << "DIVISION TEST FAILED: " << std::endl;
      std::cerr << "a:= " << polyToString(a) << ";" << std::endl;
      std::cerr << "b:= " << polyToString(b ) << ";" << std::endl;
      std::cerr << "r:= " << polyToString(r) << ";" << std::endl;
      std::cerr << "retErr " <<  retErr << std::endl;
      exit(1);
    }

    std::cerr << "Division test passed" << std::endl;
    

    
    freePolynomial_AA(ca);

  

   }

////////////////////////////////////////////////////////////////


int main(void) {


   AltArr_t* aa=makePolynomial_AA(10);
   AltArr_t* bb=makePolynomial_AA(10);

     
   mpq_t a;
   mpq_init(a);
   mpq_set_ui(a,1,5);

   ratNum_t k;
   mpq_init(k);
   mpq_set(k,a);
   addTerm_AA(aa,2 ,k);
  // addTerm_AA(aa,2,k);
   //addTermSafe_AA(aa,1,k);

  sortPolynomial_AA(aa);
  bb=deepCopyPolynomial_AA(aa);

  
 AAElem_t*c=multiplyTerms_AA(&(aa->elems[1]),&(bb->elems[1]));

  int d=c->deg;
  
   mpq_set(k,c->coef);

  addTerm_AA(bb,d,k);
  negatePolynomial_AA(bb);
  //printAA(aa);
  // printAA(bb) ;
  AltArr_t* cc= addPolynomials_AA( aa,  bb);
    //AltArr_t*cc= subPolynomials_AA( aa,bb);
   // printAA(aa);
//AltArr_t* cc=  multiplyPolynomials_AA(aa,bb);

// AltArr_t* q = NULL;
// AltArr_t* r = NULL;
//  cc=addPolynomials_AA( cc,  aa) ;
// dividePolynomials_AA(aa,cc,&q,&r);

// printAA(cc);
// printf("quotient:");
// printf("aa");

// printAA(aa);
// printf("cc");


// printAA(q);
// printf("reminder:");
// printAA(r);

// printf("\n\n\n\n\n");

// 	//int g =add(5);
// //	printf("answer is %d ", g);

time_t seed = time(NULL);

// ExpressionTree tree=convertToExpressionTree( q);


//testAdditionAndMultiplication(seed, 200, 10, 3, 1);
//testDivision(seed, 6, 10, 3, 1);

/*ExpressionTree tree;*/
// std::cout<<tree.toMapleString()<<std::endl;
// std::cout<<polyToString(q)<<std::endl;

	return 0;
}
