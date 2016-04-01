//
//  pthread_barrier.h
//  Simpson-4.0
//
//  Created by Thomas Vosegaard on 9/8/13.
//  Copyright (c) 2013 Thomas Vosegaard. All rights reserved.
//

#if defined(__APPLE__) || defined(DEFINE_PTHREAD_BARRIER_ROUTINES)
#ifndef PTHREAD_BARRIER_H_
#define PTHREAD_BARRIER_H_

#include <pthread.h>
#include <errno.h>

typedef int pthread_barrierattr_t;
typedef struct
{
    pthread_mutex_t mutex;
    pthread_cond_t cond;
    int count;
    int tripCount;
} pthread_barrier_t;


int pthread_barrier_init(pthread_barrier_t *barrier, const pthread_barrierattr_t *attr, unsigned int count);
int pthread_barrier_destroy(pthread_barrier_t *barrier);
int pthread_barrier_wait(pthread_barrier_t *barrier);


#endif // PTHREAD_BARRIER_H_
#endif // __APPLE__