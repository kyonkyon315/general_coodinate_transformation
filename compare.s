	.file	"test3.cpp"
	.text
#APP
	.globl _ZSt21ios_base_library_initv
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC8:
	.string	"=== Result ===\n"
.LC9:
	.string	"Template version : "
.LC10:
	.string	" s\n"
.LC11:
	.string	"Naive for-if ver : "
.LC12:
	.string	"Speedup ratio    : "
.LC13:
	.string	"x\n"
.LC14:
	.string	"sum1="
.LC15:
	.string	" sum2="
#NO_APP
	.section	.text.startup,"ax",@progbits
	.p2align 4
	.globl	main
	.type	main, @function
main:
.LFB2759:
	.cfi_startproc
	endbr64
	pushq	%r14
	.cfi_def_cfa_offset 16
	.cfi_offset 14, -16
	pushq	%r13
	.cfi_def_cfa_offset 24
	.cfi_offset 13, -24
	pushq	%r12
	.cfi_def_cfa_offset 32
	.cfi_offset 12, -32
	pushq	%rbp
	.cfi_def_cfa_offset 40
	.cfi_offset 6, -40
	movl	$100000, %ebp
	pushq	%rbx
	.cfi_def_cfa_offset 48
	.cfi_offset 3, -48
	subq	$48, %rsp
	.cfi_def_cfa_offset 96
	movq	$0x000000000, 32(%rsp)
	movq	$0x000000000, 40(%rsp)
	call	_ZNSt6chrono3_V212system_clock3nowEv@PLT
	movq	%rax, %r14
.L2:
	movq	$0x000000000, 8(%rsp)
	movl	$1, %ebx
	.p2align 4,,10
	.p2align 3
.L3:
	pxor	%xmm1, %xmm1
	movsd	.LC1(%rip), %xmm0
	cvtsi2sdl	%ebx, %xmm1
	addl	$1, %ebx
	mulsd	%xmm1, %xmm0
	movsd	%xmm1, 24(%rsp)
	call	sin@PLT
	movsd	24(%rsp), %xmm1
	movsd	.LC2(%rip), %xmm4
	movsd	%xmm0, 16(%rsp)
	mulsd	%xmm1, %xmm4
	movapd	%xmm4, %xmm0
	call	cos@PLT
	mulsd	16(%rsp), %xmm0
	addsd	8(%rsp), %xmm0
	movsd	%xmm0, 8(%rsp)
	cmpl	$1000, %ebx
	jne	.L3
	movsd	32(%rsp), %xmm0
	addsd	8(%rsp), %xmm0
	movl	$1, %ebx
	movq	$0x000000000, 8(%rsp)
	movsd	%xmm0, 32(%rsp)
	.p2align 4,,10
	.p2align 3
.L4:
	pxor	%xmm1, %xmm1
	movsd	.LC3(%rip), %xmm0
	cvtsi2sdl	%ebx, %xmm1
	addl	$1, %ebx
	mulsd	%xmm1, %xmm0
	movsd	%xmm1, 24(%rsp)
	call	sin@PLT
	movsd	24(%rsp), %xmm1
	movsd	.LC4(%rip), %xmm5
	movsd	%xmm0, 16(%rsp)
	mulsd	%xmm1, %xmm5
	movapd	%xmm5, %xmm0
	call	cos@PLT
	mulsd	16(%rsp), %xmm0
	addsd	8(%rsp), %xmm0
	movsd	%xmm0, 8(%rsp)
	cmpl	$1000, %ebx
	jne	.L4
	movsd	32(%rsp), %xmm0
	addsd	8(%rsp), %xmm0
	movl	$1, %ebx
	movq	$0x000000000, 8(%rsp)
	movsd	%xmm0, 32(%rsp)
	.p2align 4,,10
	.p2align 3
.L5:
	pxor	%xmm1, %xmm1
	movsd	.LC5(%rip), %xmm0
	cvtsi2sdl	%ebx, %xmm1
	addl	$1, %ebx
	mulsd	%xmm1, %xmm0
	movsd	%xmm1, 24(%rsp)
	call	sin@PLT
	movsd	24(%rsp), %xmm1
	movsd	.LC6(%rip), %xmm6
	movsd	%xmm0, 16(%rsp)
	mulsd	%xmm1, %xmm6
	movapd	%xmm6, %xmm0
	call	cos@PLT
	mulsd	16(%rsp), %xmm0
	addsd	8(%rsp), %xmm0
	movsd	%xmm0, 8(%rsp)
	cmpl	$1000, %ebx
	jne	.L5
	movsd	32(%rsp), %xmm0
	addsd	8(%rsp), %xmm0
	movsd	%xmm0, 32(%rsp)
	subl	$1, %ebp
	jne	.L2
	call	_ZNSt6chrono3_V212system_clock3nowEv@PLT
	movl	$100000, %ebp
	movq	%rax, %r12
	call	_ZNSt6chrono3_V212system_clock3nowEv@PLT
	movq	%rax, %r13
.L7:
	movq	$0x000000000, 8(%rsp)
	movl	$1, %ebx
	.p2align 4,,10
	.p2align 3
.L8:
	pxor	%xmm1, %xmm1
	movsd	.LC1(%rip), %xmm0
	cvtsi2sdl	%ebx, %xmm1
	addl	$1, %ebx
	mulsd	%xmm1, %xmm0
	movsd	%xmm1, 24(%rsp)
	call	sin@PLT
	movsd	24(%rsp), %xmm1
	movsd	.LC2(%rip), %xmm3
	movsd	%xmm0, 16(%rsp)
	mulsd	%xmm1, %xmm3
	movapd	%xmm3, %xmm0
	call	cos@PLT
	mulsd	16(%rsp), %xmm0
	addsd	8(%rsp), %xmm0
	movsd	%xmm0, 8(%rsp)
	cmpl	$1000, %ebx
	jne	.L8
	movsd	40(%rsp), %xmm0
	addsd	8(%rsp), %xmm0
	movl	$1, %ebx
	movq	$0x000000000, 8(%rsp)
	movsd	%xmm0, 40(%rsp)
	.p2align 4,,10
	.p2align 3
.L9:
	pxor	%xmm1, %xmm1
	movsd	.LC3(%rip), %xmm0
	cvtsi2sdl	%ebx, %xmm1
	addl	$1, %ebx
	mulsd	%xmm1, %xmm0
	movsd	%xmm1, 24(%rsp)
	call	sin@PLT
	movsd	24(%rsp), %xmm1
	movsd	.LC4(%rip), %xmm3
	movsd	%xmm0, 16(%rsp)
	mulsd	%xmm1, %xmm3
	movapd	%xmm3, %xmm0
	call	cos@PLT
	mulsd	16(%rsp), %xmm0
	addsd	8(%rsp), %xmm0
	movsd	%xmm0, 8(%rsp)
	cmpl	$1000, %ebx
	jne	.L9
	movsd	40(%rsp), %xmm0
	addsd	8(%rsp), %xmm0
	movl	$1, %ebx
	movq	$0x000000000, 8(%rsp)
	movsd	%xmm0, 40(%rsp)
	.p2align 4,,10
	.p2align 3
.L10:
	pxor	%xmm1, %xmm1
	movsd	.LC5(%rip), %xmm0
	cvtsi2sdl	%ebx, %xmm1
	addl	$1, %ebx
	mulsd	%xmm1, %xmm0
	movsd	%xmm1, 24(%rsp)
	call	sin@PLT
	movsd	24(%rsp), %xmm1
	movsd	.LC6(%rip), %xmm2
	movsd	%xmm0, 16(%rsp)
	mulsd	%xmm1, %xmm2
	movapd	%xmm2, %xmm0
	call	cos@PLT
	mulsd	16(%rsp), %xmm0
	addsd	8(%rsp), %xmm0
	movsd	%xmm0, 8(%rsp)
	cmpl	$1000, %ebx
	jne	.L10
	movsd	40(%rsp), %xmm0
	addsd	8(%rsp), %xmm0
	movsd	%xmm0, 40(%rsp)
	subl	$1, %ebp
	jne	.L7
	call	_ZNSt6chrono3_V212system_clock3nowEv@PLT
	subq	%r14, %r12
	pxor	%xmm0, %xmm0
	movsd	.LC7(%rip), %xmm1
	cvtsi2sdq	%r12, %xmm0
	subq	%r13, %rax
	leaq	_ZSt4cout(%rip), %rbx
	movl	$15, %edx
	movq	%rbx, %rdi
	leaq	.LC8(%rip), %rsi
	leaq	.LC10(%rip), %r12
	divsd	%xmm1, %xmm0
	movsd	%xmm0, 8(%rsp)
	pxor	%xmm0, %xmm0
	cvtsi2sdq	%rax, %xmm0
	movapd	%xmm0, %xmm7
	divsd	%xmm1, %xmm7
	movq	%xmm7, %rbp
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	movl	$19, %edx
	leaq	.LC9(%rip), %rsi
	movq	%rbx, %rdi
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	movsd	8(%rsp), %xmm0
	movq	%rbx, %rdi
	call	_ZNSo9_M_insertIdEERSoT_@PLT
	movl	$3, %edx
	movq	%r12, %rsi
	movq	%rax, %rdi
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	movl	$19, %edx
	leaq	.LC11(%rip), %rsi
	movq	%rbx, %rdi
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	movq	%rbp, %xmm0
	movq	%rbx, %rdi
	call	_ZNSo9_M_insertIdEERSoT_@PLT
	movl	$3, %edx
	movq	%r12, %rsi
	movq	%rax, %rdi
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	movl	$19, %edx
	leaq	.LC12(%rip), %rsi
	movq	%rbx, %rdi
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	movq	%rbx, %rdi
	movq	%rbp, %xmm7
	divsd	8(%rsp), %xmm7
	movapd	%xmm7, %xmm0
	call	_ZNSo9_M_insertIdEERSoT_@PLT
	movl	$2, %edx
	leaq	.LC13(%rip), %rsi
	movq	%rax, %rdi
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	movl	$5, %edx
	leaq	.LC14(%rip), %rsi
	movq	%rbx, %rdi
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	movsd	32(%rsp), %xmm0
	movq	%rbx, %rdi
	call	_ZNSo9_M_insertIdEERSoT_@PLT
	movl	$6, %edx
	leaq	.LC15(%rip), %rsi
	movq	%rax, %rdi
	movq	%rax, %rbx
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	movsd	40(%rsp), %xmm0
	movq	%rbx, %rdi
	call	_ZNSo9_M_insertIdEERSoT_@PLT
	movq	%rax, %rbp
	movq	(%rax), %rax
	movq	-24(%rax), %rax
	movq	240(%rbp,%rax), %rbx
	testq	%rbx, %rbx
	je	.L23
	cmpb	$0, 56(%rbx)
	je	.L13
	movzbl	67(%rbx), %eax
.L14:
	movq	%rbp, %rdi
	movsbl	%al, %esi
	call	_ZNSo3putEc@PLT
	movq	%rax, %rdi
	call	_ZNSo5flushEv@PLT
	addq	$48, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 48
	xorl	%eax, %eax
	popq	%rbx
	.cfi_def_cfa_offset 40
	popq	%rbp
	.cfi_def_cfa_offset 32
	popq	%r12
	.cfi_def_cfa_offset 24
	popq	%r13
	.cfi_def_cfa_offset 16
	popq	%r14
	.cfi_def_cfa_offset 8
	ret
.L13:
	.cfi_restore_state
	movq	%rbx, %rdi
	call	_ZNKSt5ctypeIcE13_M_widen_initEv@PLT
	movq	(%rbx), %rax
	movl	$10, %esi
	movq	%rbx, %rdi
	call	*48(%rax)
	jmp	.L14
.L23:
	call	_ZSt16__throw_bad_castv@PLT
	.cfi_endproc
.LFE2759:
	.size	main, .-main
	.section	.rodata.cst8,"aM",@progbits,8
	.align 8
.LC1:
	.long	-755914244
	.long	1062232653
	.align 8
.LC2:
	.long	-755914244
	.long	1063281229
	.align 8
.LC3:
	.long	-1133871366
	.long	1063818100
	.align 8
.LC4:
	.long	-755914244
	.long	1064329805
	.align 8
.LC5:
	.long	1202590843
	.long	1064598241
	.align 8
.LC6:
	.long	-1133871366
	.long	1064866676
	.align 8
.LC7:
	.long	0
	.long	1104006501
	.ident	"GCC: (Ubuntu 13.3.0-6ubuntu2~24.04) 13.3.0"
	.section	.note.GNU-stack,"",@progbits
	.section	.note.gnu.property,"a"
	.align 8
	.long	1f - 0f
	.long	4f - 1f
	.long	5
0:
	.string	"GNU"
1:
	.align 8
	.long	0xc0000002
	.long	3f - 2f
2:
	.long	0x3
3:
	.align 8
4:
