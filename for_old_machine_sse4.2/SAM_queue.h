/*
  Authors:
  Haoyu Cheng
	  chhy@mail.ustc.edu.cn
*/



#ifndef __SAMQUENUE_
#define __SAMQUENUE_

#include "Auxiliary.h"
#include <sys/queue.h>
#include "Process_Reads.h"
#include <string>

typedef struct element_queue
{
    int ele;
}element_queue;


typedef struct Queue
{
    Read  *qBase;
    int front;
    int rear;  ///rear指向队尾的下一个元素
    int output_queue_length;
}QUEUE;

typedef struct
{
    char* sVal;
    int length;
	///std::string  stringVal;
} tmp_result_string;


typedef struct SAM_Queue
{
    tmp_result_string  *qBase;
    int front;
    int rear;  ///rear指向队尾的下一个元素
    int output_queue_length;
}SAM_Queue;




struct SINGLE_END_QUEUE{
    char* single_location_result_string;
    int length;
    TAILQ_ENTRY(SINGLE_END_QUEUE) entries;
};
TAILQ_HEAD(HEADER_SINGLE_END_QUEUE, SINGLE_END_QUEUE);

void push_SINGLE_END_QUEUE(struct HEADER_SINGLE_END_QUEUE* header, int* queue_length, char* alignemnt, int alignment_length);

char* pop_all_SINGLE_END_QUEUE(struct HEADER_SINGLE_END_QUEUE* header, int* queue_length, int total_length);

void initQueue(QUEUE *pq, int length);
void enQueue(QUEUE *pq , Read* value);
void deQueue(QUEUE *pq , Read *value);

int if_empty(QUEUE *pq);
int if_full(QUEUE *pq);

inline int if_empty(QUEUE *pq)
{
	if (pq->front == pq->rear)
		return 1;
	else
		return 0;
}

inline int if_full(QUEUE *pq)
{
	if (pq->front == (pq->rear + 1) % (pq->output_queue_length))
		return 1;
	else
		return 0;
}

inline void enQueue(QUEUE *pq, Read* value)
{

	strcpy(pq->qBase[pq->rear].name, value->name);
	strcpy(pq->qBase[pq->rear].qual, value->qual);
	strcpy(pq->qBase[pq->rear].rseq, value->rseq);
	strcpy(pq->qBase[pq->rear].seq, value->seq);
	pq->qBase[pq->rear].length = value->length;
	pq->rear = (pq->rear + 1) % (pq->output_queue_length);
}

inline void deQueue(QUEUE *pq, Read *value)
{

	strcpy(value->name, pq->qBase[pq->front].name);
	strcpy(value->qual, pq->qBase[pq->front].qual);
	strcpy(value->rseq, pq->qBase[pq->front].rseq);
	strcpy(value->seq, pq->qBase[pq->front].seq);
	value->length = pq->qBase[pq->front].length;
	pq->front = (pq->front + 1) % (pq->output_queue_length);
}


#endif
