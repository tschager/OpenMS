// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Thomas Tschager $
// $Authors: Simon Rösch and Thomas Tschager $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_DENOVO_SYMDIFF_PEPTIDESOLVERSYMMETRICDIFFERENCEGENERAL_H
#define OPENMS_ANALYSIS_DENOVO_SYMDIFF_PEPTIDESOLVERSYMMETRICDIFFERENCEGENERAL_H

#include "symdiffscore_multyions.h"
#include "config.h"
#include <iostream>
#include <cmath>
#include <set>

/**
  @brief This method is the core of our symmetric difference algorithm.

  From a given set of peaks the matrix spectrum graph is computed. Then in order to find the @em rho best solutions, a DFS is executed. To do this efficiently we first compute L and R that allow us to abort the DFS search once we are on a path that will not lead to a good enough solution.

*/
void find_all_peptides(vector<peak>& peaks, char last_letter, double rho, vector<peptide>& result, config conf)
{
	scoreCacheClear();

	int k = peaks.size();

	int edge_count=0;

	vector< vector<simple_edge> > matrix_edges = vector< vector<simple_edge> >(2*k, vector<simple_edge>() );
	vector< vector<simple_edge> > matrix_in_edges = vector< vector<simple_edge> >(2*k, vector<simple_edge>() );
	vector< simple_edge > edges;

	simple_edge end_edge;
	end_edge.from=1;
	end_edge.to=1;
	end_edge.label="";
	end_edge.id=edge_count++;
	edges.push_back(end_edge);

	//Compute Edges
	for(int x=0; x<k/2; x++)
	{
		for(int m=k-1; m>x; m--)
		{
			vector<string> label = get_label(peaks[m].mz-peaks[x].mz, conf);
			if(label.size()>0)
			{
				for(int i=0; i<label.size(); i++)
				{
					string l = label[i];
					if(l.size()>1)
					{
						double missing = peaks[x].mz+get_aminoacid_mass(l[1]);
						int l=x; int r=m;
						while(l+1<r)
						{
							int next = (l+r)/2;
							if(peaks[next].mz<missing)
							{
								l=next;
							} 
							else 
							{
								r=next;
							}
						}
						if(fabs(peaks[l].mz-missing)<conf.EPS || fabs(peaks[r].mz-missing)<conf.EPS)
						{
							continue;
						}
					}

					simple_edge e;
					e.from = 2*x;
					e.to = 2*m;
					e.label = l;
					e.id = edge_count++;
					matrix_edges[e.from].push_back(e);
					matrix_in_edges[e.to].push_back(e);
					edges.push_back(e);
				}
			}
		}
	}
	
	for(int y=k-1; y>k/2; y--)
	{
		for(int m=0; m<y; m++)
		{
			vector<string> label = get_label(peaks[y].mz-peaks[m].mz, conf);
			if(label.size()>0)
			{
				for(int i=0; i<label.size(); i++)
				{
					string l = label[i];
					//Restrict last edge to end with last_letter
					if(y!=k-1 || l[0]==last_letter || (l.size()>1 && l[l.size()-2]==last_letter))
					{
						if(l.size()>1)
						{
							double missing = peaks[m].mz+get_aminoacid_mass(l[1]);
							int l=m; int r=y;
							while(l+1<r)
							{
								int next = (l+r)/2;
								if(peaks[next].mz<missing)
								{
									l=next;
								}
								else
								{
									r=next;
								}
							}
							if( fabs(peaks[l].mz-missing)<conf.EPS || fabs(peaks[r].mz-missing)<conf.EPS )
							{
								continue;
							}
						}
						simple_edge e;
						e.from = 2*(k-1-y)+1;
						e.to = 2*(k-1-m)+1;
						e.label = l;
						e.id = edge_count++;
						matrix_edges[e.from].push_back(e);
						matrix_in_edges[e.to].push_back(e);
						edges.push_back(e);
					}
				}
			}
		}
	}

	//Compute L-Values and max_score in topological order
	vector< vector<double> > matrix_l = vector< vector<double> >(k, vector<double>(edge_count, NEG_INF));
	vector< vector<double> > matrix_r = vector< vector<double> >(k, vector<double>(edge_count, NEG_INF));

	double max_score = 0;
	matrix_l[0][0] = conf.START_SCORE; 
	for(int x=0; x<k; x++)
	{
		for(int e=0; e<edge_count; e++)
		{
			if(matrix_l[x][e]==NEG_INF)
			{
				continue;
			}
			else if(edges[e].to/2+x/2==k-1)
			{
				max_score = max(max_score, matrix_l[x][e]);
				matrix_r[x][e] = 0;
			}
			else
			{
				for(int i=0; i<matrix_edges[x].size(); i++)
				{
					simple_edge out_e = matrix_edges[x][i];
					if(edges[e].to/2+out_e.to/2>k-1)
					{
						continue;
					}
					else if(out_e.to > edges[e].to)
					{
						matrix_l[edges[e].to][out_e.id] = max(matrix_l[edges[e].to][out_e.id], matrix_l[x][e]+score(peaks, out_e, edges[e], conf));
					}
					else
					{
						matrix_l[out_e.to][edges[e].id] = max(matrix_l[out_e.to][edges[e].id], matrix_l[x][e]+score(peaks, out_e, edges[e], conf));
					}
				}
			}
		}
	}

	//Compute R-Values in backward topological order
	for(int x=k-1; x>=0; x--)
	{
		for(int e=0; e<edge_count; e++)
		{
			if(matrix_r[x][e]>NEG_INF)
			{
				for(int i=0; i<matrix_in_edges[x].size(); i++)
				{
					simple_edge in_e = matrix_in_edges[x][i];
					if(in_e.from < edges[e].from)
					{
						matrix_r[edges[e].from][in_e.id] = max(matrix_r[edges[e].from][in_e.id], matrix_r[x][e]+score(peaks, edges[e], in_e, conf));
					}
					else
					{
						matrix_r[in_e.from][edges[e].id] = max(matrix_r[in_e.from][edges[e].id], matrix_r[x][e]+score(peaks, in_e, edges[e], conf));
					}
				}
			}
		}
	}
	matrix_r[0][0] = max_score;



	//find all paths with score >= rho*max_score, for this we use DFS

	vector< int > edge_counter = vector< int >(k, 0);
	int x = 0;
	int e = 0;
	vector<simple_edge> edge_stack;
	vector<double> scores;
	vector<int> e_stack;
	vector<int> x_stack;
	scores.push_back(conf.START_SCORE);
	while (true)
	{
		if (scores.back()+matrix_r[x][e]>=rho*max_score && edge_counter[x]<matrix_edges[x].size())
		{

			if (edges[e].to/2+matrix_edges[x][edge_counter[x]].to/2>k-1)
			{
				edge_counter[x]++;
				continue;
			}

			edge_stack.push_back(matrix_edges[x][edge_counter[x]]);
			edge_counter[x]++;

			scores.push_back(scores.back() + score(peaks, edge_stack.back(), edges[e], conf));
			e_stack.push_back(e);
			x_stack.push_back(x);
			if (edge_stack.back().to<edges[e].to)
			{
				x=edge_stack.back().to;
			}
			else
			{
				x = edges[e].to;
				e = edge_stack.back().id;
			}
			edge_counter[x] = 0;

			if (edges[e].to/2+x/2==k-1 && scores.back()>=rho*max_score)
			{

				vector<string> p;
				for(int i=0; i<edge_stack.size(); i++)
				{
					simple_edge e = edge_stack[i];
					if(e.from%2==0)
					{
						p.push_back(e.label);
					}
				}
				for(int i=edge_stack.size()-1; i>=0; i--)
				{
					simple_edge e = edge_stack[i];
					if(e.from%2!=0)
					{
						p.push_back(e.label);
					}
				}

				result.push_back(mp(scores.back(),p));

			}

		}
		else
		{
			if(edge_stack.size()==0)
			{
				break;//we are done
			}
			else
			{
				x = x_stack.back();
				x_stack.pop_back();
				e = e_stack.back();
				e_stack.pop_back();
				scores.pop_back();
				edge_stack.pop_back();
			}
		}

	}
	
}

#endif
