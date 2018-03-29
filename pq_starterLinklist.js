// CPCS 324 Algorithms & Data Structures 2
// Outline - Priority queue data structure
// 2017, Dr. Muhammad Al-Hashimi


// -----------------------------------------------------------------------
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

function PQueue() {
	/**#@+
	 * @description Property field - 
	 */
	// user input fields
	/** head pointer of  priority queue
	   @private */
	this.pq = new List(); // requirement: linked-list implementation
	/**#@- */

	// specify (design) methods
	/**#@+
	 * @description <em>Member method - </em> 
	 */
	/**  return true if the PQ is empty or false if not  {uses linklist package method: .isEmpty()}*/
	this.isEmpty = isEmptyPQImpl;
	/**   remove/return item with minimum priority (weight) from the PQ  */
	this.deleteMin = deleteMinPQImpl;
	/** insert an item with priority */
	this.insert = insertPQImpl;
	/**update item priority (decrease as defined in textbook)*/
	this.decrease = decreasePQImpl;
	/**#@- */
}

// -----------------------------------------------------------------------
// Priority queue node constructor (document using JSDOC comments)

function PQNode(item, key) {
	/**#@+
	 * @description Property field - 
	 */
	// user input fields
	/** to store the item in the priority queue*/
	this.item = item;
	/**key to maintain queue order*/
	this.prior = key;
	/**#@- */

	// specify (design) methods

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
function insertPQImpl(item, key) {
	var item = {
		dataitem: item,
		key: key
	};
	var insertNode = new LNode(item);
	var hptr = this.pq.first;
	// if pq is empty create node and insert it 	
	if (this.isEmpty()) {
		this.pq.insert(item);

		// else if key of the first is greater than the inserted node key, insert the node at first of pq
	} else if (hptr.item.key > item.key) {
		insertNode.next = hptr;
		this.pq.first = insertNode;
		//else insert in the right place
	} else {
		while (hptr.next != null) {
			if (hptr.next.item.key > item.key) {
				break;
			}
			hptr = hptr.next;
		}
		insertNode.next = hptr.next;
		hptr.next = insertNode;

	}
}
/**
check if the PQ is empty 
@methodOf PQueue#
@author Hatoon Mohammed
 @returns {boolean} true if empty, false if not empty
 */
function isEmptyPQImpl() {
	return this.pq.isEmpty();
}
/**
update the item priority 
@methodOf PQueue#
@author Hatoon Mohammed
 @param {integer} item 
 @param {integer} key
 */
function decreasePQImpl(item, key) {
	// 2 help pointers
	var hptr = this.pq.first;
	var hptr2;
	// travers to the end of the PQ
	while (hptr !== null) {
		//if found the node
		if (item === hptr.item.dataitem) {
			// update the node priority
			if (hptr.item.key > key) {
				hptr.item.key = key;
			}
			// if the item is at the head of PQ
			if (hptr.item.dataitem === this.pq.first.item.dataitem) {
				//delete the head and store it's info then reinsert it to the PQ
				var itemTemp = hptr.item.dataitem;
				var keyTemp = hptr.item.key;
				var head = this.deleteMin();
				this.insert(itemTemp, keyTemp);
			} else {
				var itemTemp = hptr.item.dataitem;
				var keyTemp = hptr.item.key;
				hptr2.next = hptr.next;
				this.insert(hptr.item.dataitem, hptr.item.key);
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
function deleteMinPQImpl() {
	//check if pq is not empty
	if (!this.isEmpty()) {
		//using the linklist delete first function
		var deleted = this.pq.delete_first();
		return deleted;
	}
}