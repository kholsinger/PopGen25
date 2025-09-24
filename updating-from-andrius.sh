#!/bin/bash

# Set upstream repo URL here
UPSTREAM_URL="https://github.com/adagilis/PopGen25.git"
BRANCH="main"

# Check if upstream remote exists
if ! git remote | grep -q upstream; then
  echo "Adding upstream remote..."
  git remote add upstream "$UPSTREAM_URL"
else
  echo "Upstream remote already exists."
fi

echo "Fetching upstream changes..."
git fetch upstream

echo "Checking out $BRANCH branch..."
git checkout $BRANCH

echo "Merging upstream/$BRANCH into $BRANCH..."
git merge upstream/$BRANCH

echo "Pushing changes to origin..."
git push origin $BRANCH

echo "Fork updated successfully!"
