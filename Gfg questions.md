**1]K largest elements from a big file or array**
```
#include <vector>
#include <queue>
using namespace std;
vector<int> kLargestElements(const vector<int>& arr, int k) {
    priority_queue<int, vector<int>, greater<int>> minHeap;
    for (int num : arr) {
        minHeap.push(num);
        if (minHeap.size() > k) {
            minHeap.pop();
        }
    }
    vector<int> result;
    while (!minHeap.empty()) {
        result.push_back(minHeap.top());
        minHeap.pop();
    }
    return result;
}
```

**2]Find a triplet a, b, c such that a2 = b2 + c2**
```
#include <vector>
#include <set>
using namespace std;
bool findTriplet(const vector<int>& arr) {
    set<int> squares;
    for (int num : arr) {
        squares.insert(num * num);
    }
    for (int i = 0; i < arr.size(); ++i) {
        for (int j = i + 1; j < arr.size(); ++j) {
            int a2 = arr[i] * arr[i] + arr[j] * arr[j];
            if (squares.find(a2) != squares.end()) {
                return true; // Triplet found
            }
        }
    }
    return false; // No triplet found
}
```

**3]Binary tree traversal** <br />
Left view
```
import java.util.*;
void leftView(TreeNode root) {
    if (root == null) return;
    Queue<TreeNode> queue = new LinkedList<>();
    queue.add(root);
    while (!queue.isEmpty()) {
        int n = queue.size();
        for (int i = 0; i < n; i++) {
            TreeNode node = queue.poll();
            if (i == 0) System.out.print(node.val + " ");
            if (node.left != null) queue.add(node.left);
            if (node.right != null) queue.add(node.right);
        }
    }
}
```
Right view
```
void rightView(TreeNode root) {
    if (root == null) return;
    Queue<TreeNode> queue = new LinkedList<>();
    queue.add(root);
    while (!queue.isEmpty()) {
        int n = queue.size();
        for (int i = 0; i < n; i++) {
            TreeNode node = queue.poll();
            if (i == n - 1) System.out.print(node.val + " ");
            if (node.left != null) queue.add(node.left);
            if (node.right != null) queue.add(node.right);
        }
    }
}
```
Top view
```
void topView(TreeNode root) {
    if (root == null) return;
    Map<Integer, Integer> map = new TreeMap<>();
    Queue<Pair<TreeNode, Integer>> queue = new LinkedList<>();
    queue.add(new Pair<>(root, 0));
    while (!queue.isEmpty()) {
        Pair<TreeNode, Integer> pair = queue.poll();
        TreeNode node = pair.getKey();
        int hd = pair.getValue();
        if (!map.containsKey(hd)) map.put(hd, node.val);
        if (node.left != null) queue.add(new Pair<>(node.left, hd - 1));
        if (node.right != null) queue.add(new Pair<>(node.right, hd + 1));
    }
    for (int val : map.values()) System.out.print(val + " ");
}
```
Bottom view
```
void bottomView(TreeNode root) {
    if (root == null) return;
    Map<Integer, Integer> map = new TreeMap<>();
    Queue<Pair<TreeNode, Integer>> queue = new LinkedList<>();
    queue.add(new Pair<>(root, 0));
    while (!queue.isEmpty()) {
        Pair<TreeNode, Integer> pair = queue.poll();
        TreeNode node = pair.getKey();
        int hd = pair.getValue();
        map.put(hd, node.val);
        if (node.left != null) queue.add(new Pair<>(node.left, hd - 1));
        if (node.right != null) queue.add(new Pair<>(node.right, hd + 1));
    }
    for (int val : map.values()) System.out.print(val + " ");
}
```
Maximum of a level
```
int maxOfLevel(TreeNode root, int level) {
    if (root == null) return Integer.MIN_VALUE;
    Queue<TreeNode> queue = new LinkedList<>();
    queue.add(root);
    int currentLevel = 0;
    while (!queue.isEmpty()) {
        int n = queue.size();
        int max = Integer.MIN_VALUE;
        if (currentLevel == level) {
            for (int i = 0; i < n; i++) {
                TreeNode node = queue.poll();
                max = Math.max(max, node.val);
            }
            return max;
        }
        for (int i = 0; i < n; i++) {
            TreeNode node = queue.poll();
            if (node.left != null) queue.add(node.left);
            if (node.right != null) queue.add(node.right);
        }
        currentLevel++;
    }
    return Integer.MIN_VALUE; // Level not found
}
```
Minimum of a level
```
int minOfLevel(TreeNode root, int level) {
    if (root == null) return Integer.MAX_VALUE;
    Queue<TreeNode> queue = new LinkedList<>();
    queue.add(root);
    int currentLevel = 0;
    while (!queue.isEmpty()) {
        int n = queue.size();
        int min = Integer.MAX_VALUE;
        if (currentLevel == level) {
            for (int i = 0; i < n; i++) {
                TreeNode node = queue.poll();
                min = Math.min(min, node.val);
            }
            return min;
        }
        for (int i = 0; i < n; i++) {
            TreeNode node = queue.poll();
            if (node.left != null) queue.add(node.left);
            if (node.right != null) queue.add(node.right);
        }
        currentLevel++;
    }
    return Integer.MAX_VALUE; // Level not found
}
```
Children Sum Property
```
boolean childrenSumProperty(TreeNode node) {
    if (node == null || (node.left == null && node.right == null)) return true;
    int sum = 0;
    if (node.left != null) sum += node.left.val;
    if (node.right != null) sum += node.right.val;
    return (node.val == sum) && childrenSumProperty(node.left) && childrenSumProperty(node.right);
}
```
Diameter of a binary tree
```
int diameter(TreeNode root) {
    int[] diameter = new int[1];
    height(root, diameter);
    return diameter[0];
}
int height(TreeNode node, int[] diameter) {
    if (node == null) return 0;
    int leftHeight = height(node.left, diameter);
    int rightHeight = height(node.right, diameter);
    diameter[0] = Math.max(diameter[0], leftHeight + rightHeight);
    return 1 + Math.max(leftHeight, rightHeight);
}
```
Inorder Traversal
```
void inorderTraversal(TreeNode root) {
    if (root == null) return;
    inorderTraversal(root.left);
    System.out.print(root.val + " ");
    inorderTraversal(root.right);
}
```
Preorder Traversal
```
void preorderTraversal(TreeNode root) {
    if (root == null) return;
    System.out.print(root.val + " ");
    preorderTraversal(root.left);
    preorderTraversal(root.right);
}
```
Postorder Traversal
```
void postorderTraversal(TreeNode root) {
    if (root == null) return;
    postorderTraversal(root.left);
    postorderTraversal(root.right);
    System.out.print(root.val + " ");
}
```
Level Order Traversal
```
void levelOrderTraversal(TreeNode root) {
    if (root == null) return;
    Queue<TreeNode> queue = new LinkedList<>();
    queue.add(root);
    while (!queue.isEmpty()) {
        TreeNode node = queue.poll();
        System.out.print(node.val + " ");
        if (node.left != null) queue.add(node.left);
        if (node.right != null) queue.add(node.right);
    }
}
```

**4]Convert Binary Tree to Doubly Linked List (DLL)**
```
TreeNode prev = null;
TreeNode head = null;
void convertToDLL(TreeNode root) {
    if (root == null) return;
    convertToDLL(root.left);
    if (prev == null) {
        head = root; // This node becomes the head of the DLL
    } else {
        root.left = prev;
        prev.right = root;
    }
    prev = root;
    convertToDLL(root.right);
}
```

**5]Lowest Common ancestor in a Binary Search Tree and Binary Tree.**
Binary Search Tree
```
TreeNode lowestCommonAncestorBST(TreeNode root, TreeNode p, TreeNode q) {
    if (root == null) return null;
    if (p.val < root.val && q.val < root.val) {
        return lowestCommonAncestorBST(root.left, p, q);
    } else if (p.val > root.val && q.val > root.val) {
        return lowestCommonAncestorBST(root.right, p, q);
    } else {
        return root;
    }
}
```
Binary Tree
```
TreeNode lowestCommonAncestorBinaryTree(TreeNode root, TreeNode p, TreeNode q) {
    if (root == null || root == p || root == q) return root;
    TreeNode left = lowestCommonAncestorBinaryTree(root.left, p, q);
    TreeNode right = lowestCommonAncestorBinaryTree(root.right, p, q);
    if (left != null && right != null) {
        return root;
    }
    return (left != null) ? left : right;
}
```

**6]Implement a stack with push(), pop() and min() in O(1) time.**
```
#include <stack>
#include <stdexcept>
using namespace std;
class MinStack {
    stack<int> mainStack;
    stack<int> minStack;
public:
    void push(int x) {
        mainStack.push(x);
        if (minStack.empty() || x <= minStack.top()) {
            minStack.push(x);
        }
    }
    void pop() {
        if (mainStack.empty()) throw runtime_error("Stack is empty");
        int poppedValue = mainStack.top();
        mainStack.pop();
        if (poppedValue == minStack.top()) {
            minStack.pop();
        }
    }
    int top() {
        if (mainStack.empty()) throw runtime_error("Stack is empty");
        return mainStack.top();
    }
    int min() {
        if (minStack.empty()) throw runtime_error("Stack is empty");
        return minStack.top();
    }
};
```

**7]Reverse Linked List in Groups of Size k**
```
class ListNode {
    int val;
    ListNode next;
    ListNode(int x) { val = x; }
}
class Solution {
    public ListNode reverseKGroup(ListNode head, int k) {
        ListNode current = head;
        ListNode prevTail = null;
        ListNode newHead = null;
        while (current != null) {
            ListNode groupHead = current;
            ListNode prev = null;
            int count = 0;
            while (current != null && count < k) {
                current = current.next;
                count++;
            }
            if (count == k) {
                current = groupHead;
                while (count > 0) {
                    ListNode nextNode = current.next;
                    current.next = prev;
                    prev = current;
                    current = nextNode;
                    count--;
                }
                if (newHead == null) {
                    newHead = prev; // Set the new head for the first reversed group
                }

                if (prevTail != null) {
                    prevTail.next = prev; // Connect previous tail to the new head of this group
                }
                prevTail = groupHead; // Set the previous tail for the next group
            } else {
                if (prevTail != null) {
                    prevTail.next = groupHead; // Connect previous tail to the remaining nodes
                }
                break; // No more groups to reverse
            }
        }
        return newHead != null ? newHead : head; // Return the new head or the original head
    }
}
```

**8]Given two numbers represented by two linked lists, write a function that returns sum list**
```
class ListNode {
    int val;
    ListNode next;
    ListNode(int x) { val = x; }
}
class Solution {
    public ListNode addTwoNumbers(ListNode l1, ListNode l2) {
        ListNode dummyHead = new ListNode(0);
        ListNode current = dummyHead;
        int carry = 0;
        while (l1 != null || l2 != null || carry != 0) {
            int sum = carry;
            if (l1 != null) {
                sum += l1.val;
                l1 = l1.next;
            }
            if (l2 != null) {
                sum += l2.val;
                l2 = l2.next;
            }
            carry = sum / 10;
            current.next = new ListNode(sum % 10);
            current = current.next;
        }
        return dummyHead.next;
    }
}
```

**9]Rotate a matrix by 90 degree.**
```
void rotateMatrix(vector<vector<int>>& matrix) {
    int n = matrix.size();
    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            swap(matrix[i][j], matrix[j][i]);
        }
    }
    for (int i = 0; i < n; ++i) {
        reverse(matrix[i].begin(), matrix[i].end());
    }
}
```

**10] Stock span problem**
```
vector<int> calculateSpan(const vector<int>& prices) {
    int n = prices.size();
    vector<int> span(n);
    stack<int> s;
    for (int i = 0; i < n; ++i) {
        while (!s.empty() && prices[s.top()] <= prices[i]) {
            s.pop();
        }
        span[i] = (s.empty()) ? (i + 1) : (i - s.top());
        s.push(i);
    }
    return span;
}
```

**11]Next greater element**
```
vector<int> nextGreaterElement(const vector<int>& nums) {
    int n = nums.size();
    vector<int> result(n, -1);
    stack<int> s;
    for (int i = 0; i < n; ++i) {
        while (!s.empty() && nums[s.top()] < nums[i]) {
            result[s.top()] = nums[i];
            s.pop();
        }
        s.push(i);
    }
    return result;
}
```

**12]Maximum sum subarray such that no elements are consecutive**
```
int maxSumNonConsecutive(const vector<int>& nums) {
    int n = nums.size();
    if (n == 0) return 0;
    if (n == 1) return nums[0];
    vector<int> dp(n, 0);
    dp[0] = nums[0];
    dp[1] = max(nums[0], nums[1]);
    for (int i = 2; i < n; ++i) {
        dp[i] = max(dp[i - 1], nums[i] + dp[i - 2]);
    }
    return dp[n - 1];
}
```

**13]Edit distance**
```
int editDistance(const string& word1, const string& word2) {
    int m = word1.size();
    int n = word2.size();
    vector<vector<int>> dp(m + 1, vector<int>(n + 1));
    for (int i = 0; i <= m; ++i) {
        for (int j = 0; j <= n; ++j) {
            if (i == 0) {
                dp[i][j] = j; 
            } else if (j == 0) {
                dp[i][j] = i;
            } else if (word1[i - 1] == word2[j - 1]) {
                dp[i][j] = dp[i - 1][j - 1]; 
            } else {
                dp[i][j] = 1 + min({dp[i - 1][j], dp[i][j - 1], dp[i - 1][j - 1]});
            }
        }
    }
    return dp[m][n];
}
```

**14]Assembly Scheduling**
```
int assemblyLineScheduling(const vector<vector<int>>& a, const vector<vector<int>>& t, int e1, int e2, int x1, int x2) {
    int n = a[0].size(); 
    vector<int> dp1(n), dp2(n);
    dp1[0] = e1 + a[0][0];
    dp2[0] = e2 + a[1][0];
    for (int i = 1; i < n; ++i) {
        dp1[i] = min(dp1[i - 1] + a[0][i], dp2[i - 1] + t[1][i] + a[0][i]);
        dp2[i] = min(dp2[i - 1] + a[1][i], dp1[i - 1] + t[0][i] + a[1][i]);
    }
    return min(dp1[n - 1] + x1, dp2[n - 1] + x2);
}
```
