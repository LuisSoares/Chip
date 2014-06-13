

#!/usr/bin/env python
''' Author:Luis Soares
    Version:0.1
    Create plots of chip-seq data
'''
#importing necessary modules: os,pickle,tkinter,numpy,np,time,matplotlib,collections
import os
import pickle

try:
    from tkinter import Tk,filedialog
except:
    print('Please install tkinter to run application')
try:
    import numpy as np
except:
    print('Please install numpy to run application')

import time

try:
    import matplotlib.pyplot as plt
except:
    print('Please install matplotlib to run application')

from matplotlib.patches import Rectangle
from collections import OrderedDict

#end of imports

#Setting up global variables to be used

#color scheme list to be used in graphical plotting
color_scheme=[(127.5/255, 0, 0),
              (237/255,28/255,36/255),
              (241/255,140/255,34/255),
              (255/255, 222/255, 23/255),
              (173/255, 209/255, 54/255),
              (8/255, 135/255 ,67/255),
              (71/255, 195/255, 211/255),
              (33/255, 64/255, 154/255),
              (150/255,100/255,155/255),
              (238/255, 132/255, 181/255)]

#default chromosome sizes for S.cerevisiae R64-1-1
chrom=OrderedDict([('chrI',230218),
                 ('chrII',813184),
                 ('chrIII',316620),
                 ('chrIV',1531933),
                 ('chrV',576874),
                 ('chrVI',270161),
                 ('chrVII',1090940),
                 ('chrVIII',562643),
                 ('chrIX',439888),
                 ('chrX',745751),
                 ('chrXI',666816),
                 ('chrXII',1078177),
                 ('chrXIII',924431),
                 ('chrXIV',784333),
                 ('chrXV',1091291),
                 ('chrXVI',948066)])
       
#Setting up classes to be used

class Genome:
    #class definition for what is a genome, contains 2 parameters:
    #1-Species name
    #2-Ordered dictionary containing chromosome name as key and size as value
    def __init__(self,species,chrom):
        self.species=species
        self.chrom=chrom    
    def number_of_chromosomes(self):
        #function to print number of chromosomes
        print(len(self.chrom.keys())) 
    def genome_size(self):
        #function to print length of genome
        print(sum([value for value in self.chrom.values()]))
    def print_chromosomes(self):
        #function to print chromosome name and lenght
        a=self.chrom.keys()
        for item in a:
            print(item+'\t-->\t'+str(self.chrom[item]))
        


class Chromosome:
    #class definition for what is a chromosome, contains 2 parameters:
    #1-name of the chromosome
    #2-size of the chromosome
    def __init__(self,name,size):
        self.name=name
        self.size=size

class Chromosome_track(Chromosome):
    #class definition for what is a chromosome track,
    #inherits from Chromosome class
    def __init__(self,track_name,track):
        self.track=track
        super().__init__(self.track[0],chrom[self.track[0]])
        self.starting_track=np.zeros(self.size)
        self.track=track
        self.track_name=track_name
        #changed to account for overlapping intervals
        #significantly slows down (2x) the process (even more if there are ovelapping intervals need to find better method
        #remove 'if' part(leaving the 'else' part) if sure there are no overlapping intervals
        #WARNING: THERE IS A PROBLEM if start and end in bed file are inverted
        for item in self.track[1]:            
            #if (not 0) in self.starting_track[(item[0]):(item[1])]:
             #   self.zero_track=np.zeros(self.size)
              #  self.zero_track[(item[0]):(item[1])]=item[2]
               # self.starting_track=np.add(self.starting_track,self.zero_track)
            #else:
                self.starting_track[(item[0]):(item[1])]=item[2]
    #def __repr__(self):
    #   print (self.starting_track[0:20])
    #   print (self.name)
    #   print(len(self.starting_track))
    #   return ''
   
class gene:
	def __init__(self,chromosome,code,start,end,strand,name='x'):
		self.chromosome=chromosome
		self.code=code
		self.start=int(start)
		self.end=int(end)
		self.strand=strand
		self.name=name
		self.length=self.end-self.start

# Start of funtion definitions
def get_path():
    #function to set current working path using filedialog
    while True:
        root=Tk()
        path=filedialog.askdirectory(initialdir='.')
        root.destroy()
        if os.path.exists(path):
            return path
        else:pass


#Creates a empty default genome
def get_genome():        
    print("Loading Genome")
    start=time.time()
    default_genome=Genome('S.cerevisiae',chrom)
    end=time.time()
    print ('Elapsed Time= '+str(round(end-start,2))+'sec')
    return default_genome

#loads a gene list
def load_gene_list(g_list):
    #change open to context dependent
    print ('Loading genes list')
    start=time.time()
    genes=open(g_list)
    genes.readline()
    gene_list=[]
    for line in genes:
        a=line.split('\t')
        temp1=gene(a[0],a[4].rstrip(),a[1],a[2],a[3],a[5].rstrip())
        gene_list.append(temp1)
    genes.close()
    end=time.time()
    print ('Elapsed Time= '+str(round(end-start,2))+'sec')
    return gene_list


def extract_tracks(file,chromosome,filetype='mo'):#expand for other filetypes
    '''returns a tupple containg in first position the number of the chromosome
    followed by a list of tupples(start,end,value)'''
    
    track=[]#track will be a list of tuple(start,end,value)
    with open(file) as input_data:
        input_data.readline() #discard header line
        for a in input_data:
            b=a.split('\t')
            if b[0]==chromosome:
                track.append((float(b[1]),float(b[2]),float(b[3])))
    return (chromosome,track)



genome_tracks=OrderedDict() 
def get_genome_tracks():
    root=Tk()
    input_file=filedialog.askopenfilename(initialdir=path)
    root.destroy()
    data_set_name=''
    while data_set_name=='':
        data_set_name=input('Enter name for Data set: ')
    
    print ("Loading Data")
    start=time.time()
    genome_tracks_temp=[]
    for item in chrom.keys():
            genome_tracks_temp.append(Chromosome_track(data_set_name,extract_tracks(input_file,item)))
    genome_tracks[data_set_name]=genome_tracks_temp
    end=time.time()
    print ('Elapsed Time= '+str(round(end-start,2))+'sec')

def plot_gene(gene_range,extend_range,i):
    #function getting too big, needs spliting
       
        color_number=0
        
        max_y=0
        if i==0: #added to not plot automatically (use in top4 function)
            plt.figure(figsize=(7,5))
            plt.subplots_adjust(bottom=0.12)
            plot1=plt.subplot(1,1,1)
        for chr_item in genome_tracks.keys():
            for item in genome_tracks[chr_item]:
                if item.name==gene_range[2]:
                    
                    if gene_range[3]=='+':
                        plt.xlim(int(gene_range[0])-extend_range,int(gene_range[1])+extend_range)
                        x=np.array(range(int(gene_range[0])-extend_range,int(gene_range[1])+extend_range))
                        y=item.starting_track[int(gene_range[0])-extend_range:int(gene_range[1])+extend_range]
                        plt.plot(x,y,'-',color=color_scheme[color_number],linewidth=1.5,aa=True,label=chr_item)
                        color_number+=1
                        if max(y)>max_y:max_y=max(y)
                        
                    
                    else:
                        plt.xlim(int(gene_range[1])+extend_range,int(gene_range[0])-extend_range)
                        x=np.array(range(int(gene_range[1])+extend_range,int(gene_range[0])-extend_range,-1))
                        y=item.starting_track[int(gene_range[1])+extend_range:int(gene_range[0])-extend_range:-1]
                        plt.plot(x,y,'-',color=color_scheme[color_number],linewidth=1.5,aa=True,label=chr_item)
                        color_number+=1
                        if max(y)>max_y:max_y=max(y)
                        
        if extend_range!=0:
                            plot1.add_patch(Rectangle((int(gene_range[0]), 0), int(gene_range[1])-int(gene_range[0]), max_y*1.1, fc="blue",alpha=0.1))            
        plt.legend(fontsize=14,loc='best')
        if i==0:
            plot1.set_title(gene_range[-1]+'/'+gene_range[-2],fontsize=16,fontstyle='italic')
            plot1.spines['bottom'].set_linewidth('1.5')
            plot1.spines['top'].set_linewidth('1.5')
            plot1.spines['right'].set_linewidth('1.5')
            plot1.spines['left'].set_linewidth('1.5')
            plot1.yaxis.set_tick_params(labelsize=13)
            plot1.xaxis.set_tick_params(labelsize=12)
            plot1.set_xlabel(gene_range[2],fontsize=14)
            plot1.set_ylabel('Enrichment',fontsize=14)
        plt.ylim(0,max_y*1.1)
        if i==0:
            plt.show()

def get_gene(i=0,gene=''):#arguments added so to be able to plot automatically if
    gene_range=()         #arguments are passed, i is only here to not plot if !=0
    extend_range=0
    if gene=='':
        user_input=input('Enter gene: ').upper()
    else:user_input=gene
    for item in gene_list:
        if item.name==user_input.rstrip() or item.code==user_input.rstrip():
            gene_range=(item.start,item.end,item.chromosome,item.strand,item.code,item.name)
    if not gene_range:
        print('No such gene\n')
    if i==0:
        extend_range=input('Extend range: ')
    if not extend_range:
        extend_range=0
    else:
        extend_range=int(extend_range)
    plot_gene(gene_range,extend_range,i)

def anchor_plot():
    #the anchor plot exludes cases where the window size causes the plot to fall
    #off chromosome
    color_number=0
    answer=int(input('Enter Option (1-centered range,2-skewed range): '))
    if answer==1:
        lower_range=-int(input('Enter range: '))
        higher_range=-lower_range
    elif answer==2:
        lower_range=int(input('Lower range: '))
        higher_range=int(input('Higher range: '))
    else:
        print('Invalid option')
        return
    if lower_range>higher_range:print('Unusable ranges')
    plt.figure(figsize=(6,6))
    plt.subplots_adjust(bottom=0.12)
    plot1=plt.subplot(1,1,1)
    window=abs(lower_range)+abs(higher_range)
    max_height=0
    min_height=1000000
    for __ in genome_tracks.keys():#ugly ugly ugly
        track=genome_tracks[__]
        x=np.array(range(lower_range,higher_range))
        start=time.time()
        y=np.zeros(len(range(lower_range,higher_range)))
        n=0
        for item in gene_list:
            gene_range=(item.start,item.end,item.chromosome,item.strand,item.code)
            for _ in track:
                if _.name==gene_range[2]:
                    try:
                        if gene_range[3]=='+':
                            temp=_.starting_track[int(gene_range[0])+lower_range:int(gene_range[0])+higher_range]
				#print(len(y))
				#print(len(temp))
                            
                            y=np.add(y,temp)
                        else:
                            temp=_.starting_track[int(gene_range[1])-higher_range:int(gene_range[1])-lower_range]
                            temp=temp[::-1]
				#print(len(y))
				#print(len(temp))
                            y=np.add(y,temp)
                            

                        n=n+1
                        
                    except:pass
        plt.plot(x,y/n,'-',color=color_scheme[color_number],linewidth=1.5,label=__)
        if max(y/n)>max_height:max_height=max(y/n)
        if min(y/n)<min_height:min_height=min(y/n)
        plt.legend(fontsize=12,loc='best')
        color_number+=1
    print(str(len(gene_list)-n)+" genes excluded from analysis (range outside chromosome length)")
    end=time.time()
    print('Elapsed Time= '+str(round(end-start,2))+'sec')
    #plt.xticks(range(0,window+1,int((window/4))),range(lower_range,higher_range+1,int((window/4))))
    plot1.spines['bottom'].set_linewidth('1.5')
    plot1.spines['top'].set_linewidth('1.5')
    plot1.spines['right'].set_linewidth('1.5')
    plot1.spines['left'].set_linewidth('1.5')
    plot1.yaxis.set_tick_params(labelsize=12)
    plot1.xaxis.set_tick_params(labelsize=12)
    plot1.set_ylabel('Average Enrichment (n=%s)'%n,fontsize=14)
    plot1.set_xlabel('0=ATG',fontsize=14)
    plt.xlim(lower_range,higher_range)
    plt.ylim(min_height*0.9,max_height*1.1)
    plt.show()

    
def heat_map(window):
    #can use this function to get a tuple containing the heatmap matrix and the number of lines
    #can be useful if other types of sorting are required 
    start=time.time()
    y=np.zeros(window*2,dtype=np.float16)
    n=0
    for item in gene_list:
        gene_range=(item.start,item.end,item.chromosome,item.strand,item.code)
        for _ in genome_tracks[list(genome_tracks.keys())[0]]:
            if _.name==gene_range[2]:
                try:
                    if gene_range[3]=='+':
                        temp=_.starting_track[int(gene_range[0])-window:int(gene_range[0])+window]
				#print(len(y))
				#print(len(temp))
                        y=np.vstack((y,temp))
                    else:
                        temp=_.starting_track[int(gene_range[1])-window:int(gene_range[1])+window]
                        temp=temp[::-1]
                        y=np.vstack((y,temp))
                    
                    n=n+1
                except:pass
    end=time.time()
    print ('Elapsed Time= '+str(round(end-start,2))+'sec')
    return y,n
    

def plot_heat_map():
    Range=int(input('Enter range: '))
    values,n=heat_map(Range)
    plt.figure(figsize=(5,7))
    plt.subplots_adjust(left=0.17)
    plot1=plt.subplot(1,1,1)
    plt.ylim(0,len(values))
    plt.xticks(range(0,2*Range+1,int((Range*2/4))),range(-Range,Range+1,int((Range*2/4))))
    plot1.spines['bottom'].set_linewidth('1.5')
    plot1.spines['top'].set_linewidth('1.5')
    plot1.spines['right'].set_linewidth('1.5')
    plot1.spines['left'].set_linewidth('1.5')
    plot1.yaxis.set_tick_params(labelsize=12)
    plot1.xaxis.set_tick_params(labelsize=12)
    values=values[np.argsort(np.array([np.mean(item) for item in values]))] #sorts by average                                      #to sort
    plt.pcolormesh(values,cmap=plt.cm.Reds,rasterized=True)#rasterize makes it faster to manipulate but less precise display
    plt.colorbar(shrink=0.25,aspect=5)
    plot1.set_ylabel('Genes (n=%s)'%n,fontsize=14)
    plot1.set_xlabel('ATG',fontsize=14)
    del values
    plt.show()
    
def list_dataset():
    print('\nData sets:')
    for n,name in enumerate(genome_tracks.keys(),start=1):
        print ('\t {0}- {1}'.format(n,name))
    print()
    
def remove_dataset():
    while True:
        genome_tracks.popitem()
    

def get_metagene():
    color_number=0
    plt.figure(figsize=(6,6))
    plt.subplots_adjust(bottom=0.12)
    plot1=plt.subplot(1,1,1)
    max_height=0
    min_height=1000000
    for __ in genome_tracks.keys():#ugly ugly ugly
        track=genome_tracks[__]
        x=np.array(range(0,46))
        start=time.time()
        y=np.zeros(46)
        number_of_genes=0
        error=0
        bins=46
        for item in gene_list:
            mean=0
            for _ in track:
                if _.name==item.chromosome:
                    try:
                        if item.strand=='+':
                            gene_track=_.starting_track[item.start:item.end]
                            try:
                                slices = np.linspace(0, len(gene_track), bins+1, True).astype(np.int)
                                counts = np.diff(slices)
                                mean = np.add.reduceat(gene_track, slices[:-1]) / counts
                                #print(mean[0:10])

                                y=y+mean
                                #print(y[0:10])
                            except Exception as detail:
                                print(detail)
                                error+=1
                            number_of_genes+=1
                        else:
                            gene_track=_.starting_track[item.end:item.start:-1]
                            try:
                                slices = np.linspace(0, len(gene_track), bins+1, True).astype(np.int)
                                counts = np.diff(slices)
                                mean = np.add.reduceat(gene_track, slices[:-1]) / counts
                                #print(mean[0:10])

                                y=y+mean
                                #print(y[0:10])
                            except Exception as detail:
                                print(detail)
                                error+=1
                            number_of_genes+=1
                    except:pass
        plt.plot(x,y/number_of_genes,'-',color=color_scheme[color_number],linewidth=1.5,label=__)
        plt.legend(fontsize=12,loc='best')
        color_number+=1
    print(str(len(gene_list)-number_of_genes)+" genes excluded from analysis (range outside chromosome length)")
    end=time.time()
    print('Elapsed Time= '+str(round(end-start,2))+'sec')
    plot1.spines['bottom'].set_linewidth('1.5')
    plot1.spines['top'].set_linewidth('1.5')
    plot1.spines['right'].set_linewidth('1.5')
    plot1.spines['left'].set_linewidth('1.5')
    plot1.yaxis.set_tick_params(labelsize=12)
    plot1.xaxis.set_tick_params(labelsize=12)
    plt.xticks([0,11.5,22,33.5,44.999], [0,25,50,75,100])
    plot1.set_ylabel('Average Enrichment (n=%s)'%number_of_genes,fontsize=14)
    plot1.set_xlabel('0=ATG',fontsize=14)
    plt.xlim(0,45)
    plt.show()

def save_tracks():
    #saves loaded tracks as a pickle file containing just the genome_tracks
    #dictionary
    with open('tracks.pkl','wb') as output:
        start=time.time()
        pickle.dump(genome_tracks,output)
        end=time.time()
        print('Elapsed Time= '+str(round(end-start,2))+'sec')

def load_compiled_tracks():
    #if a file named tracks.pkl is present it replaces genome tracks with the
    #compiled trackds
    try:
        with open('tracks.pkl','rb') as input_file:
            start=time.time()
            loaded=pickle.load(input_file)
            end=time.time()
            print('Elapsed Time= '+str(round(end-start,2))+'sec')
    except:
        print('Compiled tracks file not found')
    return loaded

def plot_chromosome():
    a=[np.reshape(genome_tracks[list(genome_tracks.keys())[0]][n].starting_track,
    (len(genome_tracks[list(genome_tracks.keys())[0]][n].starting_track),1)) for n in range(0,
    len(genome_tracks[list(genome_tracks.keys())[0]]))]
    left, width = 0.06, 0.025
    bottom=0.1
    plt.figure(figsize=(10,5))
    temp=[]
    for position,n in enumerate(a):
	    left_h = left+0.055*(position)+0.02
	    chromosome= [left_h, bottom, width, (len(n)/1531933)*0.8]
	    temp.append(plt.axes(chromosome))
	    plt.ylim(0,len(n))
	    temp[position].axes.get_yaxis().set_visible(False)
	    temp[position].axes.get_xaxis().set_visible(False)
	    temp[position].pcolormesh(n,cmap=plt.cm.Reds,rasterized=True)
	    temp[position].spines['bottom'].set_linewidth('2')
	    temp[position].spines['top'].set_linewidth('2')
	    temp[position].spines['right'].set_linewidth('2')
	    temp[position].spines['left'].set_linewidth('2')
	    temp[position].axes.set_title(str(position+1))
    temp[0].text(12,1500000,'Chromosome distribution',fontsize=16)
    plt.show()

def interface():
    while True:
        print('Options: ')
        list_of_options=['1-Load dataset','2-Plot specific gene','3-Plot anchor plot'
                         ,'4-Plot heat map (restricted to first Dataset)'
                         ,'5-Plot metagene analysis','6-Plot chromosome distribution (restricted to first Dataset)'
                         ,'7-List datasets','8-Remove datasets','9-Save datasets'
                         ,'0-Load compiled tracks','Q-Quit']
        for item in list_of_options:
            print('\t{}'.format(item))
        user_input=input('Option: ')
        try:
            if user_input=='1':
                get_genome_tracks()
            elif user_input=='2':
                get_gene()
            elif user_input=='3':
                anchor_plot()
            elif user_input=='4':
                plot_heat_map()
            elif user_input=='5':
                get_metagene()
            elif user_input=='6':
                plot_chromosome()
            elif user_input=='7':
                list_dataset()
            elif user_input=='8':
                remove_dataset()
            elif user_input=='9':
                save_tracks()
            elif user_input=='0':
                global genome_tracks
                genome_tracks=load_compiled_tracks()
            elif user_input.upper()=='Q':
                return
            else:
                print('Invalid option')
        except:
            print('Error occured')


        
#Start of default loading functions        
#establishing current working path where gene list file should be found
path=get_path()        
os.chdir(path)
default_genome=get_genome()
gene_list=load_gene_list('genes_tss.mo')    

#if __name__ == "__main__":
interface()

def plot_region():
    chromosome=input('Chromosome: ')
    start=input('Start: ')
    end=input('End: ')
    for chr_item in genome_tracks.keys():
        for item in genome_tracks[chr_item]:

