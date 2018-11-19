/*
 * grad.h
 *
 * Written 2014-2015 by Tuomo Valkonen <tuomov@iki.fi>, University of Cambridge.
 *
 */

#include <string.h>
//#include <libtu/types.h>
#include "number.h"

#ifndef _GRAD_H
#define _GRAD_H

static __inline__ FNumber gradx(size_t h, size_t w, const FNumber (*u)[h][w],
                                size_t y, size_t x)
{
    // Neumann boundary, forward differences
    return (x<w-1 ? (*u)[y][x+1]-(*u)[y][x] : 0);
}

static __inline__ FNumber grady(size_t h, size_t w, const FNumber (*u)[h][w],
                                size_t y, size_t x)
{
    // Neumann boundary, forward differences
    return (y<h-1 ? (*u)[y+1][x]-(*u)[y][x] : 0);
}

static __inline__ FNumber gradx_transpose(size_t h, size_t w, const FNumber (*u)[h][w],
                                          size_t y, size_t x)
{
    // Neumann boundary, forward differences
    return (x==0 ? -(*u)[y][x] 
                 : (x==w-1 ? (*u)[y][x-1]
                           : (*u)[y][x-1]-(*u)[y][x]));
}

static __inline__ FNumber grady_transpose(size_t h, size_t w, const FNumber (*u)[h][w],
                                          size_t y, size_t x)
{
    // Neumann boundary, forward differences
    return (y==0 ? -(*u)[y][x] 
                 : (y==h-1 ? (*u)[y-1][x]
                           : (*u)[y-1][x]-(*u)[y][x]));
}

static __inline__ void divupdatex(size_t h, size_t w, FNumber (*u)[h][w], 
                                  size_t y, size_t x, FNumber v)
{
    // Neumann boundary, forward differences
    if(x<w-1){
        (*u)[y][x]-=v;
        (*u)[y][x+1]+=v;
    }
}

static __inline__ void divupdatey(size_t h, size_t w, FNumber (*u)[h][w], 
                                  size_t y, size_t x, FNumber v)
{
    // Neumann boundary, forward differences
    if(y<h-1){
        (*u)[y][x]-=v;
        (*u)[y+1][x]+=v;
    }
}

static __inline__ FNumber divergence(size_t h, size_t w, const FNumber (*phi)[2][h][w], 
                                     size_t y, size_t x)
{
    return -(gradx_transpose(h, w, &((*phi)[0]), y, x)
             +grady_transpose(h, w, &((*phi)[1]), y, x));
}

static __inline__ void do_grad_fn(size_t h, size_t w, 
                                  const FNumber (*u)[h][w],
                                  FNumber (*v)[2][h][w])
{
    #pragma omp parallel for
    for(size_t y=0; y<h; y++){
        for(size_t x=0; x<w; x++){
            (*v)[0][y][x]=gradx(h, w, u, y, x);
            (*v)[1][y][x]=grady(h, w, u, y, x);
        }
    }
}

static __inline__ void do_gradconj_fn(size_t h, size_t w, 
                                      const FNumber (*u)[2][h][w],
                                      FNumber (*v)[h][w])
{
    #pragma omp parallel for
    for(size_t y=0; y<h; y++){
        for(size_t x=0; x<w; x++){
            (*v)[y][x]= gradx_transpose(h, w, &((*u)[0]), y, x)
                       +grady_transpose(h, w, &((*u)[1]), y, x);
        }
    }
}

static __inline__ FNumber ax(size_t h, size_t w, const FNumber (*u)[h][w],
                             size_t y, size_t x)
{
    // counterpart to gradx
    return (x<w-1 ? (*u)[y][x+1]+(*u)[y][x] : 0);
}

static __inline__ FNumber ay(size_t h, size_t w, const FNumber (*u)[h][w],
                             size_t y, size_t x)
{
    // counterpart to grady
    return (y<h-1 ? (*u)[y+1][x]+(*u)[y][x] : 0);
}


static __inline__ void do_a_fn(size_t h, size_t w, 
                               const FNumber (*u)[h][w],
                               FNumber (*v)[h][w])
{
    #pragma omp parallel for
    for(size_t y=0; y<h; y++){
        for(size_t x=0; x<w; x++){
            (*v)[y][x]=4*(ax(h, w, u, y, x)+ay(h, w, u, y, x));
        }
    }
}

#endif /* _GRAD_H */
