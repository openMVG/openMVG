/*
	Multi-precision real number class. C++ interface fo MPFR library.
	Project homepage: http://www.holoborodko.com/pavel/
	Contact e-mail:   pavel@holoborodko.com

	Copyright (c) 2008-2011 Pavel Holoborodko

	Core Developers: 
	Pavel Holoborodko, Dmitriy Gubanov, Konstantin Holoborodko. 

	Contributors:
	Brian Gladman, Helmut Jarausch, Fokko Beekhof, Ulrich Mutze, 
	Heinz van Saanen, Pere Constans, Peter van Hoof, Gael Guennebaud, 
	Tsai Chia Cheng, Alexei Zubanov.

	****************************************************************************
	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	This library is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
	Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public
	License along with this library; if not, write to the Free Software
	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

	****************************************************************************
	****************************************************************************
	Redistribution and use in source and binary forms, with or without
	modification, are permitted provided that the following conditions
	are met:

	1. Redistributions of source code must retain the above copyright
	notice, this list of conditions and the following disclaimer.

	2. Redistributions in binary form must reproduce the above copyright
	notice, this list of conditions and the following disclaimer in the
	documentation and/or other materials provided with the distribution.

	3. The name of the author may be used to endorse or promote products
	derived from this software without specific prior written permission.

	THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
	ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
	IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
	ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
	FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
	DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
	OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
	HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
	LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
	OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
	SUCH DAMAGE.
*/
#include <cstring>
#include "mpreal.h"

#if defined (MPREAL_HAVE_CUSTOM_MPFR_MALLOC)
#include "dlmalloc.h"
#endif

using std::ws;
using std::cerr;
using std::endl;
using std::string;
using std::ostream;
using std::istream;

namespace mpfr{

mp_rnd_t   mpreal::default_rnd  = MPFR_RNDN;	//(mpfr_get_default_rounding_mode)();	
mp_prec_t  mpreal::default_prec = 64;			//(mpfr_get_default_prec)();	
int		   mpreal::default_base = 10;
int        mpreal::double_bits = -1;

#if defined (MPREAL_HAVE_CUSTOM_MPFR_MALLOC)
bool       mpreal::is_custom_malloc = false;
#endif

// Default constructor: creates mp number and initializes it to 0.
mpreal::mpreal() 
{ 

#if defined (MPREAL_HAVE_CUSTOM_MPFR_MALLOC)
	set_custom_malloc();
#endif

	mpfr_init2(mp,default_prec); 
	mpfr_set_ui(mp,0,default_rnd);

	MPREAL_MSVC_DEBUGVIEW_CODE;
}

mpreal::mpreal(const mpreal& u) 
{

#if defined (MPREAL_HAVE_CUSTOM_MPFR_MALLOC)
	set_custom_malloc();
#endif

	mpfr_init2(mp,mpfr_get_prec(u.mp));
	mpfr_set(mp,u.mp,default_rnd);

	MPREAL_MSVC_DEBUGVIEW_CODE;
}

mpreal::mpreal(const mpfr_t u)
{

#if defined (MPREAL_HAVE_CUSTOM_MPFR_MALLOC)
	set_custom_malloc();
#endif

	mpfr_init2(mp,mpfr_get_prec(u));
	mpfr_set(mp,u,default_rnd);
	
	MPREAL_MSVC_DEBUGVIEW_CODE;
}

mpreal::mpreal(const mpf_t u)
{

#if defined (MPREAL_HAVE_CUSTOM_MPFR_MALLOC)
	set_custom_malloc();
#endif

	mpfr_init2(mp,(mp_prec_t) mpf_get_prec(u)); // (gmp: mp_bitcnt_t) unsigned long -> long (mpfr: mp_prec_t)
	mpfr_set_f(mp,u,default_rnd);

	MPREAL_MSVC_DEBUGVIEW_CODE;
}

mpreal::mpreal(const mpz_t u, mp_prec_t prec, mp_rnd_t mode)
{

#if defined (MPREAL_HAVE_CUSTOM_MPFR_MALLOC)
	set_custom_malloc();
#endif

	mpfr_init2(mp,prec);
	mpfr_set_z(mp,u,mode);
	
	MPREAL_MSVC_DEBUGVIEW_CODE;
}

mpreal::mpreal(const mpq_t u, mp_prec_t prec, mp_rnd_t mode)
{

#if defined (MPREAL_HAVE_CUSTOM_MPFR_MALLOC)
	set_custom_malloc();
#endif

	mpfr_init2(mp,prec);
	mpfr_set_q(mp,u,mode);
	
	MPREAL_MSVC_DEBUGVIEW_CODE;
}

mpreal::mpreal(const double u, mp_prec_t prec, mp_rnd_t mode)
{

#if defined (MPREAL_HAVE_CUSTOM_MPFR_MALLOC)
	set_custom_malloc();
#endif

    if(double_bits == -1 || fits_in_bits(u, double_bits))
    {
    	mpfr_init2(mp,prec);
	    mpfr_set_d(mp,u,mode);
		
		MPREAL_MSVC_DEBUGVIEW_CODE;
    }
    else
        throw conversion_overflow();
}

mpreal::mpreal(const long double u, mp_prec_t prec, mp_rnd_t mode)
{ 

#if defined (MPREAL_HAVE_CUSTOM_MPFR_MALLOC)
	set_custom_malloc();
#endif

    mpfr_init2(mp,prec);
	mpfr_set_ld(mp,u,mode);
	
	MPREAL_MSVC_DEBUGVIEW_CODE;
}

mpreal::mpreal(const unsigned long int u, mp_prec_t prec, mp_rnd_t mode)
{ 

#if defined (MPREAL_HAVE_CUSTOM_MPFR_MALLOC)
	set_custom_malloc();
#endif

	mpfr_init2(mp,prec);
	mpfr_set_ui(mp,u,mode);
	
	MPREAL_MSVC_DEBUGVIEW_CODE;
}

mpreal::mpreal(const unsigned int u, mp_prec_t prec, mp_rnd_t mode)
{ 

#if defined (MPREAL_HAVE_CUSTOM_MPFR_MALLOC)
	set_custom_malloc();
#endif

	mpfr_init2(mp,prec);
	mpfr_set_ui(mp,u,mode);

	MPREAL_MSVC_DEBUGVIEW_CODE;
}

mpreal::mpreal(const long int u, mp_prec_t prec, mp_rnd_t mode)
{ 

#if defined (MPREAL_HAVE_CUSTOM_MPFR_MALLOC)
	set_custom_malloc();
#endif

	mpfr_init2(mp,prec);
	mpfr_set_si(mp,u,mode);

	MPREAL_MSVC_DEBUGVIEW_CODE;
}

mpreal::mpreal(const int u, mp_prec_t prec, mp_rnd_t mode)
{ 

#if defined (MPREAL_HAVE_CUSTOM_MPFR_MALLOC)
	set_custom_malloc();
#endif

	mpfr_init2(mp,prec);
	mpfr_set_si(mp,u,mode);

	MPREAL_MSVC_DEBUGVIEW_CODE;
}

#if defined (MPREAL_HAVE_INT64_SUPPORT)
mpreal::mpreal(const uint64_t u, mp_prec_t prec, mp_rnd_t mode)
{

#if defined (MPREAL_HAVE_CUSTOM_MPFR_MALLOC)
	set_custom_malloc();
#endif

	mpfr_init2(mp,prec);
	mpfr_set_uj(mp, u, mode); 

	MPREAL_MSVC_DEBUGVIEW_CODE;
}

mpreal::mpreal(const int64_t u, mp_prec_t prec, mp_rnd_t mode)
{

#if defined (MPREAL_HAVE_CUSTOM_MPFR_MALLOC)
	set_custom_malloc();
#endif

	mpfr_init2(mp,prec);
	mpfr_set_sj(mp, u, mode); 
	
	MPREAL_MSVC_DEBUGVIEW_CODE;
}
#endif

mpreal::mpreal(const char* s, mp_prec_t prec, int base, mp_rnd_t mode)
{

#if defined (MPREAL_HAVE_CUSTOM_MPFR_MALLOC)
	set_custom_malloc();
#endif

	mpfr_init2(mp,prec);
	mpfr_set_str(mp, s, base, mode); 

	MPREAL_MSVC_DEBUGVIEW_CODE;
}

mpreal::mpreal(const std::string& s, mp_prec_t prec, int base, mp_rnd_t mode)
{

#if defined (MPREAL_HAVE_CUSTOM_MPFR_MALLOC)
	set_custom_malloc();
#endif

	mpfr_init2(mp,prec);
	mpfr_set_str(mp, s.c_str(), base, mode); 

	MPREAL_MSVC_DEBUGVIEW_CODE;
}

mpreal::~mpreal() 
{ 
	mpfr_clear(mp);
}                           

// Operators - Assignment
mpreal& mpreal::operator=(const char* s)
{
	mpfr_t t;

#if defined (MPREAL_HAVE_CUSTOM_MPFR_MALLOC)
	set_custom_malloc();
#endif

	if(0==mpfr_init_set_str(t,s,default_base,default_rnd))
	{
		// We will rewrite mp anyway, so flash it and resize
		mpfr_set_prec(mp,mpfr_get_prec(t)); 
		mpfr_set(mp,t,mpreal::default_rnd);
		mpfr_clear(t);

		MPREAL_MSVC_DEBUGVIEW_CODE;

	}else{
		mpfr_clear(t);
	}

	return *this;
}

const mpreal fma (const mpreal& v1, const mpreal& v2, const mpreal& v3, mp_rnd_t rnd_mode)
{
	mpreal a;
	mp_prec_t p1, p2, p3;

	p1 = v1.get_prec(); 
	p2 = v2.get_prec(); 
	p3 = v3.get_prec(); 

	a.set_prec(p3>p2?(p3>p1?p3:p1):(p2>p1?p2:p1));

	mpfr_fma(a.mp,v1.mp,v2.mp,v3.mp,rnd_mode);
	return a;
}

const mpreal fms (const mpreal& v1, const mpreal& v2, const mpreal& v3, mp_rnd_t rnd_mode)
{
	mpreal a;
	mp_prec_t p1, p2, p3;

	p1 = v1.get_prec(); 
	p2 = v2.get_prec(); 
	p3 = v3.get_prec(); 

	a.set_prec(p3>p2?(p3>p1?p3:p1):(p2>p1?p2:p1));

	mpfr_fms(a.mp,v1.mp,v2.mp,v3.mp,rnd_mode);
	return a;
}

const mpreal agm (const mpreal& v1, const mpreal& v2, mp_rnd_t rnd_mode)
{
	mpreal a;
	mp_prec_t p1, p2;

	p1 = v1.get_prec(); 
	p2 = v2.get_prec(); 

	a.set_prec(p1>p2?p1:p2);

	mpfr_agm(a.mp, v1.mp, v2.mp, rnd_mode);

	return a;
}

const mpreal sum (const mpreal tab[], unsigned long int n, mp_rnd_t rnd_mode)
{
	mpreal x;
	mpfr_ptr* t;
	unsigned long int i;

	t = new mpfr_ptr[n];
	for (i=0;i<n;i++) t[i] = (mpfr_ptr)tab[i].mp;
	mpfr_sum(x.mp,t,n,rnd_mode);
	delete[] t;
	return x;
}

const mpreal remquo (long* q, const mpreal& x, const mpreal& y, mp_rnd_t rnd_mode)
{
	mpreal a;
	mp_prec_t yp, xp;

	yp = y.get_prec(); 
	xp = x.get_prec(); 

	a.set_prec(yp>xp?yp:xp);

	mpfr_remquo(a.mp,q, x.mp, y.mp, rnd_mode);

	return a;
}

template <class T>
std::string toString(T t, std::ios_base & (*f)(std::ios_base&))
{
	std::ostringstream oss;
	oss << f << t;
	return oss.str();
}

#if (MPFR_VERSION >= MPFR_VERSION_NUM(2,4,0))

std::string mpreal::toString(const std::string& format) const
{
	char *s = NULL;
	string out;

	if( !format.empty() )
	{
		if(!(mpfr_asprintf(&s,format.c_str(),mp) < 0))
		{
			out = std::string(s);

			mpfr_free_str(s);
		}
	}

	return out;
}

#endif

std::string mpreal::toString(int n, int b, mp_rnd_t mode) const
{
  (void)b;
  (void)mode;
#if (MPFR_VERSION >= MPFR_VERSION_NUM(2,4,0))

	// Use MPFR native function for output
	char format[128];
	int digits;

	digits = n > 0 ? n : bits2digits(mpfr_get_prec(mp));

	sprintf(format,"%%.%dRNg",digits);		// Default format

	return toString(std::string(format));

#else

	char *s, *ns = NULL; 
	size_t slen, nslen;
	mp_exp_t exp;
	string out;

#if defined (MPREAL_HAVE_CUSTOM_MPFR_MALLOC)
	set_custom_malloc();
#endif

	if(mpfr_inf_p(mp))
	{ 
		if(mpfr_sgn(mp)>0) return "+Inf";
		else			   return "-Inf";
	}

	if(mpfr_zero_p(mp)) return "0";
	if(mpfr_nan_p(mp))  return "NaN";

	s  = mpfr_get_str(NULL,&exp,b,0,mp,mode);
	ns = mpfr_get_str(NULL,&exp,b,n,mp,mode);

	if(s!=NULL && ns!=NULL)
	{
		slen  = strlen(s);
		nslen = strlen(ns);
		if(nslen<=slen) 
		{
			mpfr_free_str(s);
			s = ns;
			slen = nslen;
		}
		else {
			mpfr_free_str(ns);
		}

		// Make human eye-friendly formatting if possible
		if (exp>0 && static_cast<size_t>(exp)<slen)
		{
			if(s[0]=='-')
			{
				// Remove zeros starting from right end
				char* ptr = s+slen-1;
				while (*ptr=='0' && ptr>s+exp) ptr--; 

				if(ptr==s+exp) out = string(s,exp+1);
				else		   out = string(s,exp+1)+'.'+string(s+exp+1,ptr-(s+exp+1)+1);

				//out = string(s,exp+1)+'.'+string(s+exp+1);
			}
			else
			{
				// Remove zeros starting from right end
				char* ptr = s+slen-1;
				while (*ptr=='0' && ptr>s+exp-1) ptr--; 

				if(ptr==s+exp-1) out = string(s,exp);
				else		     out = string(s,exp)+'.'+string(s+exp,ptr-(s+exp)+1);

				//out = string(s,exp)+'.'+string(s+exp);
			}

		}else{ // exp<0 || exp>slen
			if(s[0]=='-')
			{
				// Remove zeros starting from right end
				char* ptr = s+slen-1;
				while (*ptr=='0' && ptr>s+1) ptr--; 

				if(ptr==s+1) out = string(s,2);
				else		 out = string(s,2)+'.'+string(s+2,ptr-(s+2)+1);

				//out = string(s,2)+'.'+string(s+2);
			}
			else
			{
				// Remove zeros starting from right end
				char* ptr = s+slen-1;
				while (*ptr=='0' && ptr>s) ptr--; 

				if(ptr==s) out = string(s,1);
				else	   out = string(s,1)+'.'+string(s+1,ptr-(s+1)+1);

				//out = string(s,1)+'.'+string(s+1);
			}

			// Make final string
			if(--exp)
			{
				if(exp>0) out += "e+"+mpfr::toString<mp_exp_t>(exp,std::dec);
				else 	  out += "e"+mpfr::toString<mp_exp_t>(exp,std::dec);
			}
		}

		mpfr_free_str(s);
		return out;
	}else{
		return "conversion error!";
	}
#endif
}


//////////////////////////////////////////////////////////////////////////
// I/O
ostream& operator<<(ostream& os, const mpreal& v)
{
	return os<<v.toString(static_cast<int>(os.precision()));
}

istream& operator>>(istream &is, mpreal& v)
{
	string tmp;
	is >> tmp;
	mpfr_set_str(v.mp, tmp.c_str(),mpreal::default_base,mpreal::default_rnd);
	return is;
}


#if defined (MPREAL_HAVE_CUSTOM_MPFR_MALLOC)
	// Optimized dynamic memory allocation/(re-)deallocation.
	void * mpreal::mpreal_allocate(size_t alloc_size)
	{
		return(dlmalloc(alloc_size));
	}

	void * mpreal::mpreal_reallocate(void *ptr, size_t old_size, size_t new_size)
	{
		return(dlrealloc(ptr,new_size));
	}

	void mpreal::mpreal_free(void *ptr, size_t size)
	{
		dlfree(ptr);
	}

	inline void mpreal::set_custom_malloc(void)
	{
		if(!is_custom_malloc)
		{
			mp_set_memory_functions(mpreal_allocate,mpreal_reallocate,mpreal_free);
			is_custom_malloc = true;
		}
	}
#endif

}

