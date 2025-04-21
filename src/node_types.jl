#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2025 Jimmy Dubuisson <jimmy.dubuisson@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it should be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#

# Binary tree node type declarations
abstract type AbstractNode end

struct EmptyNodeType <: AbstractNode end
const EmptyNode = EmptyNodeType()

mutable struct Node{T} <: AbstractNode
    key::T
    left::AbstractNode
    right::AbstractNode
end 