# This file is automatically generated from sp2graph setup.py
# Git information (specific commit, etc.)
git_revision = '0434304ad39d4101d5c8761f1a485dfacd8b59d8'
git_revision_short = git_revision[:7]
git_count = 0

# Version information
major   = 0
minor   = 0
micro   = 1

# Release tag
release = 'v'+'.'.join(map(str,[major, minor, micro]))

# Version (release + count)
version = release
if git_count > 0:
    # Add git-revision to the version string
    version += '+' + str(git_count)

# Extensive version description
label   = '0.0.1'
