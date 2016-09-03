/*                  sample.c                    */

/* Ewens sampling formula, freq are non-decreasing ,or age ordered 
and null terminated.
*/
double Psample(int *freq,double theta,int ageflag) {
    int i = 0,nsample = 0,count = 1;
    double prob = 1.0;
    if(ageflag) {
        for(i = 0;freq[i];i++) ;
        while(--i) {
            nsample += freq[i];
            prob *= theta/nsample;
        }
    nsample += freq[0];    
    } 
    else {      
        if(!freq[0]) return 0.0;
        nsample = freq[0];
        prob = 1.0/(double)freq[0];
        while(freq[++i]) {
            if(freq[i] != freq[i-1]) {
                count = 0;
                prob *= theta;
            }
        count++;
        nsample += freq[i];
        prob /= (double)(freq[i]*count);
        }
    }
    for(i = 1;i < nsample;i++)
        prob /= (1.0 + theta/(double)i);   
    if(ageflag) return prob;    
    else return prob*(double)nsample;    
}
