#!/usr/local/share/miniconda3/bin/python
import argparse,ffmpeg

def create_video(string,framerate=4,output='out.mp4'):
    """
    Usage: create_video(string,framerate=4,output='out.mp4')
    """
    output_options={
    'crf':20,
    'vf':"scale=trunc(iw/2)*2:trunc(ih/2)*2",
    'vcodec':'libx264',
    'pix_fmt':'yuv420p',
    }
    (
    ffmpeg.input(string,pattern_type='glob',framerate=framerate).output(output,**output_options).run(overwrite_output=True)
    )

def main():
    framerate= 4
    output = 'out.mp4'
    parser = argparse.ArgumentParser()
    parser.add_argument("-r","--framerate",type=int,help="Frame rate. Default is 4.")
    parser.add_argument("-i","--input",type=str,help="Wildcard of images enclose in quotations. Must end with type",required=True)
    parser.add_argument("-o","--output",type=str,help="Path and name of file. Default is out.mp4")
    args = parser.parse_args()
    if args.framerate is not None:
        framerate=args.framerate
    if args.output is not None:
        output=args.output
    input = args.input
    parser.parse_args()
    create_video(input,framerate,output)

if __name__=="__main__":
    main()
