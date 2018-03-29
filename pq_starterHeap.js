// CPCS 324 Algorithms & Data Structures 2
// Outline - Priority queue data structure
// 2017, Dr. Muhammad Al-Hashimi
// -----------------------------------------------------------------------
// Raghad's comments
// Basic design decisions and implementation planning (objects & interfaces):
// Initial Requirements:
// - Use the same interface in previous PQ, so it supports Prim 2
//
// Design Decisions:
// - based on a Heap structure.
//
//
// Using PQ functions based on a heap:
// - insert(): inserting a node in the PQ using Hinsert() function of Heap.
// - deleteMin(): remove and return the node with the least priority 
// 			      (it'll always be the root node, Using deleteRoot() function of Heap).
// - isEmpty(): check if the PQ is empty (use the HisEmpty of the Heap).
// 
// Unique functions:
// - makeMinimum(boolean: ): To decide whether to use heap as a min or max heap
//						     (Using makeMin(boolean: ) function of Heap).	
//
// Basic design decisions and implementation planning (objects & interfaces)
// initial requirements: to quickly support Dijkstra and second Prim's algorithms, 
// implement minimum priority functionality
// design decisions:
// based on the 324 linked list implementation
// how will the PQ ops be implemented?
//we chose the 4th choice, we'll manipulating linked list internals from the PQ package.
// some added PriorityQueue functions to use is defines bellow 
//insert(): inserting a node in the right place in PQ (Manipulating linklist insert function).
//deleteMin(): remove and return the node with the least priority (it'll always be the first node, we'll use the delete_first() function of linklist).
//decrease(): update the priority of a node if it's changed.
//isEmpty(): check if the PQ is empty (use the isEmpty of the linklist).
// code plan: start to write your API specifications (JSDOC comments) and think about 
// design consequences (design impact)
// Impact analysis:
//1. violating the encapsulation princeple in the insert function of the linklist, 
//because we'll be directly manipulating it from the PQ package.
//2. reusing the linklist isEmpty to check whether the PQ is empty or not.
//3. reusing the linklist delete_first to delete the first node which is also always the minimmum.
// -----------------------------------------------------------------------
// Priority queue object constructor (document using JSDOC comments)
function HPQueue()
{
    /**#@+
     * @description Property field - 
     */
    // user input fields
    /** head pointer of  priority queue
       @private */
    this.hpq = new Heap(); // requirement: Heap implementation
    /**#@- */

    //	this.minHeap = true;
    // specify (design) methods
    /**#@+
     * @description <em>Member method - </em> 
     */
    /**  return true if the PQ is empty or false if not  {uses linklist package method: .isEmpty()}*/
    this.hpisEmpty = isEmptyHPQImpl;
    /**   remove/return item with minimum priority (weight) from the PQ  */
    this.hpdeleteMin = deleteMinHPQImpl;
    /** insert an item with priority */
    this.hpinsert = insertHPQImpl;
    /**update item priority (decrease as defined in textbook)*/
    this.hpdecrease = decreaseHPQImpl;
    /**#@- */
    this.makeMinimum = makeMinimumImpl;

}

// -----------------------------------------------------------------------
// Priority queue node constructor (document using JSDOC comments)

function HPQNode(item, key)
{
    /**#@+
     * @description Property field - 
     */
    // user input fields
    /** to store the item in the priority queue*/
    this.item = item;
    /**key to maintain queue order*/
    this.prior = key;
    /**#@- */

    // specify (design) methodsPQNode

}

// -----------------------------------------------------------------------
// functions used by PQueue() object methods
// specify interface information (JSDOC comments)
// ....
/**
insert an item to the PQ in a decreasing order
@methodOf PQueue#
@author Hatoon Mohammed
 @param {integer} item to insert
 @param {integer} key the priority of the item
 */
function insertHPQImpl(item, key)
{
    return this.hpq.Hinsert(key, item);
}
/**
check if the PQ is empty 
@methodOf PQueue#
@author Hatoon Mohammed
 @returns {boolean} true if empty, false if not empty
 */
function isEmptyHPQImpl()
{
    return this.hpq.HisEmpty();
}
/**
update the item priority 
@methodOf PQueue#
@author Hatoon Mohammed
 @param {integer} item 
 @param {integer} key
 */
function decreaseHPQImpl(key, item)
{
    // 2 help pointers
    var hptr = this.hpq.first;
    var hptr2;
    // travers to the end of the PQ
    while (hptr !== null)
    {
        //if found the node
        if (item === hptr.item.dataitem)
        {
            // update the node priority
            if (hptr.item.key > key)
            {
                hptr.item.key = key;
            }
            // if the item is at the head of PQ
            if (hptr.item.dataitem === this.hpq.first.item.dataitem)
            {
                //delete the head and store it's info then reinsert it to the PQ
                var itemTemp = hptr.item.dataitem;
                var keyTemp = hptr.item.key;
                var head = this.hpdeleteMin();
                this.hpinsert(keyTemp, itemTemp);
            }
            else
            {
                var itemTemp = hptr.item.dataitem;
                var keyTemp = hptr.item.key;
                hptr2.next = hptr.next;
                this.hpinsert(hptr.item.key, hptr.item.dataitem);
            }
            break;
        }
        hptr2 = hptr;
        hptr = hptr.next;
    }
}
/**
delete the item with the highest priority, 
the list is sorted due to the PQ structure so delete first. 
@methodOf PQueue#
@author Hatoon Mohammed
 @returns item with the minimum weight
 */
function deleteMinHPQImpl()
{
    return this.hpq.deleteRoot();
}

/**
Set behaviour of Heap, true for min heap, and false for max heap
@methodOf PQueue#
@author Raghad Salem
@param {boolean} b
 */
function makeMinimumImpl(b)
{
    this.hpq.makeMin(b);
}


