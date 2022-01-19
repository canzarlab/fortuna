#ifndef RBT_CPP
#define RBT_CPP

#include "RBT.h"

void RBT::insert(string id, int s, int e) 
{
    insert(Node(id, s, e, nil));
}

void RBT::insert(Node n)
{
	n.p = nil;

    Node *y = nil, *x = root;
    Node *z = new Node(n);
    z->l = z->r = nil;

    while (x != nil)
        y = x, x = *z < *x ? x->l : x->r;
    z->p = y;
    if (y == nil)
        root = z;
    else
        (*z < *y ? y->l : y->r) = z;
    
    z->color = 1;
    insertRB(z);
} 

Node* RBT::search(Node* ptr, int x) 
{
    if (ptr == nil)
        return nil;
    if (ptr->contains(x))
        return ptr;
    if (x < ptr->s)
        return search(ptr->l, x);
    return search(ptr->r, x); 
}

void RBT::remove(Node* z) 
{ 
    if (z == nil) return; 
    Node* y = (z->l == nil || z->r == nil) ? z : successor(z); 
    Node* x = y->l != nil ? y->l : y->r;
    x->p = y->p; 
    if (y->p == nil)
        root = x;
    else 
       (y == y->p->l ? y->p->l : y->p->r) = x;
    if (z != y) (*z) = (*y);  
    if (!(y->color)) removeRB(x); 
    //if (y != nil) delete y;     
}

Node* RBT::min(Node* ptr) 
{
	if (ptr == nil) return nil;
    while (ptr->l != nil)
        ptr = ptr->l;
    return ptr;
}

Node* RBT::max(Node* ptr) 
{
	if (ptr == nil) return nil;
    while (ptr->r != nil)
        ptr = ptr->r;
    return ptr;
}

Node* RBT::successor(Node* ptr) 
{
    if (ptr->r != nil) return min(ptr->r);
    Node* y = ptr->p;
    while (y != nil && ptr == y->r)
        ptr = y, y = y->p;
    return y;
}

Node* RBT::predecessor(Node* ptr) 
{
    if (ptr->l != nil) return max(ptr->l);
    Node* y = ptr->p;
    while (y != nil && ptr == y->l)
        ptr = y, y = y->p;
    return y;
}

void RBT::clear(Node* ptr) 
{
    if (ptr == nil) return;
    clear (ptr->l);
    clear (ptr->r);
    delete ptr;  
}

void RBT::inorder(Node* ptr)
{
    if (ptr == nil) return;
    inorder(ptr->l);
    cout << ptr->id << " " << ptr->s << " " << ptr->e << endl;
    inorder(ptr->r);
}

void RBT::rotateL(Node* x) 
{
    Node* y = x->r;
    x->r = y->l;
    if (y->l != nil)
        y->l->p = x;
    y->p = x->p;
    if (x->p == nil) 
        root = y;
    else
       (x == x->p->l ? x->p->l : x->p->r) = y;
    y->l = x;
    x->p = y;
}

void RBT::rotateR(Node* y) 
{
    Node* x = y->l;    
    y->l = x->r;
    if (x->r != nil) 
        x->r->p = y;
    x->p = y->p;
    if (y->p == nil)
        root = x;
    else
        (y == y->p->l ? y->p->l : y->p->r) = x;
    x->r = y;
    y->p = x;   
}

void RBT::insertRB(Node* z)
{  
    while (z->p->color)
    {
        bool b = z->p == z->p->p->l;       
        Node* y = b ? z->p->p->r : z->p->p->l;
        if (y->color)
        {
            z->p->color = 0;
            y->color = 0;
            z->p->p->color = 1;
            z = z->p->p;
        }
        else 
        {
            if (z == (b ? z->p->r : z->p->l)) 
            {
                z = z->p;
                b ? rotateL(z) : rotateR(z);
            }
            z->p->color = 0;
            z->p->p->color = 1;
            b ? rotateR(z->p->p) : rotateL(z->p->p);
        }
    }
    root->color = 0;
}

void RBT::removeRB(Node* x) 
{ 
    while (x != root && !(x->color))
    {
        bool b = x == x->p->l;         
        Node* w = b ? x->p->r : x->p->l; 
        if (w->color)
        {
            w->color = 0;
            x->p->color = 1;
            b ? rotateL(x->p) : rotateR(x->p);
            w = (b ? x->p->r : x->p->l);
        }
        if (!(w->l->color) && !(w->r->color))
        {
            w->color = 1; 
            x = x->p;
        }
        else if (!(b ? w->r->color : w->l->color))
        {     
            (b ? w->l->color : w->r->color) = 0;
            //w->l->color = 0;
            w->color = 1;
            b ? rotateR(w) : rotateL(w);
            w = b ? x->p->r : x->p->l;
        }
        else
        {  
            w->color = x->p->color;
            x->p->color = 0;
            (b ? w->r->color : w->l->color) = 0;
            b ? rotateL(x->p) : rotateR(x->p);
            x = root;
        }            
    }
    x->color = 0; 
}

#endif
