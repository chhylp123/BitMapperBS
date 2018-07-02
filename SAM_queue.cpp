/*
  Authors:
  Haoyu Cheng
	  chhy@mail.ustc.edu.cn
*/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include<stdint.h>
#include <xmmintrin.h>
#include <mmintrin.h>
#include "SAM_queue.h"


///header是队列头
///queue_length是队列长度，这个需要更新
///alignemnt单个匹配位置的那个结果字符串；alignemnt必须是malloc出来的字符串，不能是预先定义长度的char数组
///alignment_length是alignemnt的字符串长度
void push_SINGLE_END_QUEUE(struct HEADER_SINGLE_END_QUEUE* header, int* queue_length, char* alignemnt, int alignment_length)
{



    (*queue_length)++;

    struct SINGLE_END_QUEUE *buffer_result;

	buffer_result = (SINGLE_END_QUEUE*)malloc(sizeof(struct SINGLE_END_QUEUE));

    ///alignemnt不应该包括'\0';alignment_length也不包括
    buffer_result->single_location_result_string=(char*)malloc(sizeof(char)*(alignment_length));

    memcpy(buffer_result->single_location_result_string,alignemnt,alignment_length);


    buffer_result->length=alignment_length;

    TAILQ_INSERT_TAIL(header,buffer_result,entries);

}

///header是队列头
///queue_length是队列长度，这个顺手给他更新到0吧
///total_length应该是字符串的总长度
///返回值的char*就是这一个reads的所有的结果的字符串，以\0结尾；这个是在这个函数里面malloc的空间;
///作为返回值的这个字符串也是必须要释放掉的
char* pop_all_SINGLE_END_QUEUE(struct HEADER_SINGLE_END_QUEUE* header, int* queue_length, int total_length)
{
    struct SINGLE_END_QUEUE *buffer_result;

    char* total_string=(char*)malloc(sizeof(char)*(total_length+1));

    int current_length=0;

    while((*queue_length)>0)
    {
        buffer_result=header->tqh_first;
        TAILQ_REMOVE(header,buffer_result,entries);
        (*queue_length)--;



        memcpy(total_string+current_length,buffer_result->single_location_result_string,buffer_result->length);
        current_length=current_length+buffer_result->length;
        free(buffer_result->single_location_result_string);
        free(buffer_result);
    }
    total_string[current_length]='\0';
    return total_string;
}

int SAM_Queue_if_empty(SAM_Queue *pq)
{
    if(pq->front == pq->rear)
        return 1;
    else
        return 0;
}

int SAM_Queue_if_full(SAM_Queue *pq)
{
    if(pq->front ==  (pq->rear+1)%(pq->output_queue_length))
        return 1;
    else
        return 0;
}


void SAM_Queue_initQueue(SAM_Queue *pq, int length)
{
    pq->output_queue_length=length+1;
	pq->qBase = (tmp_result_string*)malloc(sizeof(tmp_result_string)*(pq->output_queue_length));
    pq->front = pq->rear = 0;
}


void SAM_Queue_enQueue(SAM_Queue *pq , char* value, int Used_length)
{
    pq->qBase[pq->rear].sVal=value;
    pq->qBase[pq->rear].length=Used_length;
    pq->rear = (pq->rear+1)%(pq->output_queue_length);
}

void SAM_Queue_deQueue(SAM_Queue *pq , tmp_result_string *value)
{
    value->length=pq->qBase[pq->front].length;
    value->sVal=pq->qBase[pq->front].sVal;
    pq->front=(pq->front+1)%(pq->output_queue_length);
}






void initQueue(QUEUE *pq, int length)
{
    pq->output_queue_length=length+1;
	pq->qBase = (Read*)malloc(sizeof(Read)*(pq->output_queue_length));
    pq->front = pq->rear = 0;
}



