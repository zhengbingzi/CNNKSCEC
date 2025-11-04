#!/bin/bash
resolutions=("5000")
input_files=(
           
            "/mnt/sdf/zhengbingzi/CNNKSCEC/results/predictions_chr20-22_cnnklsccecabm+p0.2t4l0.003_student12.bedpe"

 )

thres=("0.97")
# ѭ������ÿ������
for res in "${resolutions[@]}"; do
    for thre in "${thres[@]}"; do
        # ��ʼ��һ�����������洢����ļ�·��
        output_files=()
        for input_file in "${input_files[@]}"; do
            # ��ȡ�ļ�����������·������ɾ����չ��
            file_name="${input_file##*/}"
            file_name_no_extension="${file_name%.*}"
            echo $file_name_no_extension
            # ��������ļ�·��
            output_file="${input_file%/*}/Peakachucluster_${file_name_no_extension}_loop${thre}.bedpe"
    
            # ��ӡ���ڴ������ļ���Ϣ
            #echo "Processing $input_file at resolution $res..."
    
            # ��������
            peakachu pool -r "$res" -i "$input_file" -o "$output_file" -t ${thre}
    
            #echo "Output saved to $output_file"

            # ������ļ�·�����ӵ�������
            output_files+=("$output_file")
        done

        # �ϲ���������ļ���һ���ļ�
        output_dir="${input_files[0]%/*}"
        merged_output_file="${output_dir}/Peakchucluster_chrsome_5kb_loop${thre}.bedpe"
        cat "${output_files[@]}" > "$merged_output_file"
        echo "Merged output saved to $merged_output_file"
        # ����ϲ��ļ�������
        num_lines=$(wc -l < "$merged_output_file")
        echo "Number of lines in merged output file ${thre}: $num_lines"
        # �������ϲ��ļ��������������ļ�����
        renamed_output_file="${output_dir}/new_sum${num_lines}_chrsome_loop${thre}.bedpe"
        mv "$merged_output_file" "$renamed_output_file"
        # ɾ���Ѿ��ϲ�������ļ�
        for output_file in "${output_files[@]}"; do
            rm "$output_file"
            echo "Deleted $output_file"
        done
    done
done

echo "All processing completed."

