/* 
 * File:   TPZTask.h
 * Author: Thiago
 *
 * Created on 19 de Junho de 2018, 14:08
 */

#ifndef TPZTASK_H
#define TPZTASK_H

#include <future>
#include <tpzautopointer.h>

class TPZTaskGroup;

// Helper class for representing tasks in a ThreadPool instance

class TPZTask {
public:

    TPZTask(const int priority, TPZAutoPointer<std::packaged_task<void(void)>> &task, TPZTaskGroup *taskGroup = NULL);
    int priority() const;
    virtual void start();
    virtual void Cancel();
    virtual ~TPZTask();

    friend class TPZTaskOrdering;
    friend class TPZThreadPool;
    
protected :
    enum EProcessingState {
        CREATED,
        SCHEDULED,
        STARTED,
        FINISHED
    };
    
    TPZAutoPointer<std::packaged_task<void(void)>> mTask;
    EProcessingState mState;
    
    TPZTaskGroup* mTaskGroup;
private:
    bool mSystemTask;
    int mPriority;
};

// Simple struct needed by std::priority_queue for ordering the items

struct TPZTaskOrdering {

    bool operator()(const TPZAutoPointer<TPZTask> &lhs, const TPZAutoPointer<TPZTask> &rhs) {
        if (lhs->mSystemTask) {
            if (rhs->mSystemTask) {
                return lhs->priority() < rhs->priority();
            } else {
                return false;
            }
        } else if (rhs->mSystemTask){
                return true;
        } else {
                return lhs->priority() < rhs->priority();
        }
    }
};


#endif /* TPZTASK_H */

