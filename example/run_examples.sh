#!/bin/bash
set -e

echo "=== Running SpaceTreeLib Examples ==="
echo ""

echo "1. Running kd_tree_example..."
./build/kd_tree_example
echo "✓ kd_tree_example completed"
echo ""

echo "2. Running orth_tree_example..."
./build/orth_tree_example
echo "✓ orth_tree_example completed"
echo ""

echo "3. Running p_tree_example..."
./build/p_tree_example
echo "✓ p_tree_example completed"
echo ""

echo "=== All examples completed successfully! ==="
